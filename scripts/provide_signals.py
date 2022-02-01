from numpy.core.numeric import base_repr
from numpy.testing._private.utils import nulp_diff
import pandas as pd
import matplotlib.pyplot as plt
import generate_colored_noises as gcn
import surogates as srg
import pandas as pd
import numpy as np
from numpy import random as rn
from scipy import signal
from scipy import io
import random as rd
import calendar
from dateutil.parser import parse

    
def noise_signal(length_n):
    '''Create noise signals'''
    return rn.uniform(0, 1, length_n)

def noise_phase():
    '''Create noise phase'''
    return rn.uniform(0, 2*np.pi)

def pink_noise_signal(length_n):
    '''Create pink noise signals'''
    return gcn.pink(length_n, state=None)

class signal_properties:
    '''
    Initialize the parameters needed to create sinthetic signals

    :param amplitude: multiplicative factor that determines the normalization of the signal, default is 1
    :param length_t: length of the temporal series needed to create the signal
    :param noise_amp: amplitude of the noise to be added
    :param freq0: numerical value of the base frequency
    :param freq1: numerical value of the other frequency, usually about one order of magnitude away from freq0 

    '''
    def __init__(self, amplitude = 1, length_t = 200, noise_amp = 0.1, freq0 = 10, freq1 = 100 ):
        self.amplitude = amplitude
        self.length_t = length_t
        self.noise_amp = noise_amp
        self.freq0 = freq0
        self.freq1 = freq1

def rossler_phase(name):

    "./data/imput/rossler_phase_r12e27N03.dat"

    sig = np.genfromtxt(name ,delimiter=',' ,
           dtype="f8", usecols=[0])
    t = np.arange(len(sig))
    return t, sig

def get_ave_values(xvalues, yvalues, n = 5):
    
    signal_length = len(xvalues)
    if signal_length % n == 0:
        padding_length = 0
    else:
        padding_length = n - signal_length//n % n
    xarr = np.array(xvalues)
    yarr = np.array(yvalues)
    xarr.resize(signal_length//n, n)
    yarr.resize(signal_length//n, n)
    xarr_reshaped = xarr.reshape((-1,n))
    yarr_reshaped = yarr.reshape((-1,n))
    x_ave = xarr_reshaped[:,0]
    y_ave = np.nanmean(yarr_reshaped, axis=1)
    
    return x_ave, y_ave

def online_ENSO_34(dataset):

    df_ENSO = pd.read_table(dataset)
    N = df_ENSO.shape[0]
    t0=1871
    dt=0.25
    time = np.arange(0, N) * dt + t0
    signal = df_ENSO.values.squeeze()

    return   time, signal


def read_ENSO_rain_manuel_files(name):
    ''' read formated file'''

    data = pd.read_csv(name).T.values
    sig = data[1]
    dates = data[2]   
    N = sig.shape[0]
    
    t0 = parse(dates[0], fuzzy=True).year
    dt=1.0#[month]/12
    time = np.arange(0, N) * dt + t0 
    #print(t0, 'dates, from:', dates[0], 'to', dates[-1], 'in intervals of',dt,' \n' )
    #t = np.arange(1871, 1871+len(sig)//12, 1/12. ) 

    return  time, sig/np.mean(sig)-1
    



def create_signal(frequencies, amplitudes, t, noise = False, gauss = False):
    ''' 
    return  a signal as the sum of cosines of different frequencies and corresponding amplitudes provided
    
    :param *kwargs:
        noise: if true noise is added, default is False
        gauss: if true a gaussian in the central period is added, default is False
        
    '''
    
    sig_cos = []
    sig_gauss =  10*np.real(np.exp(-0.001*(t-len(t)/2)**2))#*np.exp(1j*2*np.pi*2*(t-0.4)))
    sig_noise = rn.normal(0,0.25, size=(len(t)))
    sig = np.zeros(len(t))

    for i, f in enumerate(frequencies):
        sig_cos.append( np.cos( 2*np.pi * f * t) )  #+ signal.gausspulse(t - 0.4, fc=2)
        sig += amplitudes[i]*sig_cos[i]
    
    if noise: sig += sig_noise
    if gauss: sig += sig_gauss
    
    return t, sig



def multiplicative_coupling(par_sig = signal_properties):
    '''
    create coupled signal by multiplying two time series with a delay
    '''

    phase_noise = noise_phase()
    t = np.arange(0, 1, 1./par_sig.length_t)
    end = len(t)
    noise0 = noise_signal(end)
    noise1 = noise_signal(end)
    sig0 =  np.cos(np.pi*par_sig.freq0 *t + phase_noise) + noise0
    sig1 =  np.cos(np.pi*par_sig.freq1 *t + phase_noise)

    target = par_sig.amp * sig0[par_sig.delay:] * sig1[:end-par_sig.delay] + noise1
    return sig0, target


def create_victor_signal(num_segments,  base_length, delay = 2, amp_min = 10,  factor = 0.04, norm = 0.08, round_int = 3):

    '''
    Create one signal with CFC between 10 & 70 Hz
    :param sig_length: number of units of the signal
    :param directionality: imposes a directionality in the link:
       directionality = 0; No directionality
       directionality = 1; Directionality from 10Hz to 70Hz
       directionality = 2; Directionality from 70Hz to 10Hz
    :param sig_length:  signal length
    :param base_length: 8 to 12 times the maximum length of a segment, depends on norm and factor
    :param duration: duration of each segment
    '''
    
    dur = np.round(rn.rand(num_segments)*factor + norm, round_int) #duration of each temporal segment of the signal 
    freq_sequence = 1./dur #frequency of each sequence
    amp_sequence = rn.uniform(-1, 1, num_segments) + amp_min #amplitude of each sequence 
    print('dur, amplitudes sequence', dur, amp_sequence)

    π = np.pi

    def phase_signal():
        '''Create the phase signal [~10 Hz]'''
        sig_ph=[]
        t_freq=[]
        for i in range(num_segments):
            t_aux = np.arange( 0, 1.0, 1.0/int( np.round(base_length*dur[i], 2) ) )
            signal_segment = amp_sequence[i]*(np.sin(2*π*freq_sequence[i]*t_aux + 1.5*π) +1.0 )#
            t_freq.append(t_aux + i)
            sig_ph.append(signal_segment)

        signal_array = np.array([item for sublist in sig_ph for item in sublist])
        time_array = np.array([item for sublist in t_freq for item in sublist])

        plt.plot(time_array, signal_array)
        return time_array, signal_array 

    def amplitude_signal(f_amp = 70, amplitude = 10, c = 6):    
        '''Create the amplitude signal [~70 Hz]'''

        c = 6
        t = np.arange( 0, t_freq[-1], t_freq[-1]/len(t_freq) ) 
        sig_amp = []

        for i in range(len(t_freq)):
            aux1 = 1 - (1 / (1 + np.exp(-amplitude*(sig_ph[i] - c) )) )  
            aux2 = (np.sin(2*π*f_amp*t[i]) + 1)
            sig_amp.append( aux1 * aux2 ) 

        #plt.plot(t, sig_amp)
        return np.array( sig_amp )



    def directionality_signal(delay, white = True, pink = True):
        '''
        Impose directionality, delay one signal X units respect to the other
        
        :param *kwargs:
            delay: delay in system units
            white: sum white noise to the signal
            pink: sum pink noise to the signal
        '''

        #no coupled sum of the signals
        end = len(t_freq) - 1 
        no_coupling_sig = sig_ph + sig_amp

        #signals summed with delay, amp to phase
        time_delay_amp_to_phas = t_freq[delay:end]
        sig_ph_del = sig_ph[delay:end]
        sig_amp_del = sig_amp[:end-delay]
        amp_to_phase_sig = sig_ph_del + sig_amp_del

        print('end', end, 'delay', delay)

        #signals summed with delay, phase to amp
        time_delay_phas_to_amp = t_freq[delay:end]
        sig_ph_del2 = sig_ph[:end-delay]
        sig_amp_del2 = sig_amp[delay:end]
        phas_to_amp_sig = sig_ph_del2 + sig_amp_del2

        if white: 
            no_coupling_sig = no_coupling_sig + sig_white_n
            amp_to_phase_sig =  amp_to_phase_sig + sig_white_n[delay:end]
            phas_to_amp_sig =  phas_to_amp_sig + sig_white_n[:end-delay+1]

        if pink:
            no_coupling_sig = no_coupling_sig + sig_pink_n 
            amp_to_phase_sig = amp_to_phase_sig + sig_pink_n[delay:end]
            phas_to_amp_sig = phas_to_amp_sig + sig_pink_n[:end-delay+1]

        return no_coupling_sig , amp_to_phase_sig, phas_to_amp_sig, time_delay_amp_to_phas, time_delay_phas_to_amp
        
    t_freq, sig_ph = phase_signal()
    sig_amp = amplitude_signal()
    sig_white_n = noise_signal(len(t_freq))
    sig_pink_n = pink_noise_signal(len(t_freq))
    no_coupling_sig , amp_to_phase_sig, phas_to_amp_sig, time_delay_amp_to_phas, time_delay_phas_to_amp = directionality_signal(delay, white = False, pink = False)
    return t_freq, no_coupling_sig , amp_to_phase_sig, phas_to_amp_sig, time_delay_amp_to_phas, time_delay_phas_to_amp

'''compute sinusoidal signal'''
#frequencies = [ 1/20., 1/100, 1/6] #freq, they shall be below one
#amplitudes = [0.5, 1, 2]
#t = np.arange(600) 
#sampling_dt = 1

#sig = create_signal(frequencies, amplitudes, t, gauss = False )#


def plot_delayed_undelayed():

    '''ploting the cuplend and uncoupled sinthetic signals '''

    fig, ax = plt.subplots(4,1)
    ax[0].set_title("no coupling")
    ax[0].plot(t_freq, no_coupling_sig, color ='green')
    ax[1].set_title("amp to phase")
    ax[1].plot(time_delay_amp_to_phas, amp_to_phase_sig, c = 'r')
    ax[2].set_title("Phase to amp")
    ax[2].plot(time_delay_phas_to_amp, phas_to_amp_sig, c = 'b')

    s = int(len(amp_to_phase_sig)/4)
    e = int(len(amp_to_phase_sig)/3)
    print ('start, end', s, e, 'length no, am to ph, ph to am', len(no_coupling_sig), len(amp_to_phase_sig), len(phas_to_amp_sig))
    ax[3].set_title("comparison of the signals")
    ax[3].plot(t_freq[s:e], no_coupling_sig[s:e], c = 'g')
    ax[3].plot(t_freq[s:e], amp_to_phase_sig[s:e], c = 'r')
    ax[3].plot(t_freq[s:e], phas_to_amp_sig[s:e], c = 'b')
    plt.tight_layout()






#victor_sig = io.loadmat('./data/exp_pro/matlab_victor_sin_data/signal_alpha2gamma.mat')


'''compute delay signal
num_segments = 3#100
base_length = 11111
delay = 50

t_freq, no_coupling_sig , amp_to_phase_sig, phas_to_amp_sig, time_delay_amp_to_phas, time_delay_phas_to_amp = create_victor_signal(num_segments,  base_length, delay)
plot_delayed_undelayed()
plot_victor_sig_file()
plt.show()
'''



