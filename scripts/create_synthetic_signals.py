from numpy.core.numeric import base_repr
from numpy.testing._private.utils import nulp_diff
import pandas as pd
import matplotlib.pyplot as plt
import generate_colored_noises as gcn
import pandas as pd
import numpy as np
from numpy import random as rn
from scipy import signal


def create_signal(frequencies, amplitudes, t, noise = False, gauss = False):
    ''' 
    return  a signal as the sum of cosines of different frequencies and corresponding amplitudes provided
    
    :param *kwargs:
        noise: if true noise is added, default is False
        gauss: if true a gaussian in the central period is added, default is False
        
    '''
    
    sig_cos = []
    sig_gauss =  10*np.real(np.exp(-0.001*(t-len(t)/2)**2))#*np.exp(1j*2*np.pi*2*(t-0.4)))
    sig_noise = rn.normal(0,0.1, size=(len(t)))
    sig = np.zeros(len(t))

    for i, f in enumerate(frequencies):
        sig_cos.append( np.cos( 2*np.pi * f * t) )  #+ signal.gausspulse(t - 0.4, fc=2)
        sig += amplitudes[i]*sig_cos[i]
    
    if noise: sig += sig_noise
    if gauss: sig += sig_gauss
    
    return sig


def create_victor_signal(num_segments,  base_length, amp_min = 10,  factor = 0.04, norm = 0.08, round_int = 3):

    '''
    Create one signal with CFC between 10 & 70 Hz
    :param sig_length: number of units of the signal
    :param directionality: imposes a directionality in the link:
       directionality = 0; No directionality
       directionality = 1; Directionality from 10Hz to 70Hz
       directionality = 2; Directionality from 70Hz to 10Hz
    :param sig_length:  signal length
    :param Fs = 1000
    :param duration: duration of each segment
    '''
    
    dur = np.round(rn.rand(num_segments)*factor + norm, round_int) #duration of each temporal segment of the signal 
    freq_sequence = 1./dur #frequency of each sequence
    amp_sequence = rn.uniform(-1, 1, num_segments) + amp_min #amplitude of each sequence 
    print(dur, amp_sequence)

    π = np.pi

    def phase_signal():
        '''Create the phase signal [~10 Hz]'''
        sig_ph=[]
        
        for i in range(num_segments):

            print
            t_aux = np.arange(0, int( np.round(base_length*dur[i], 2) ), dur[i] ) 
            
            signal_segment = amp_sequence[i]*(np.sin(2*π*freq_sequence[i]*t_aux + 1.5*π) + 1 )
            sig_ph.append(signal_segment)

        return np.array(sig_ph)

    def amplitude_signal(f_amp = 70, amplitude = 10, c = 6):    
        '''Create the amplitude signal [~70 Hz]'''

        f_amp = 70
        a = 10
        c = 6
        t = np.arange( 0, len(sig_ph), len(sig_ph)/base_length ) 
        sig_amp = []

        for i in range(len(sig_ph)):
            aux1 = 1 - (1 / (1 + np.exp(-a*(sig_ph[i] - c) )) )  
            aux2 = (np.sin(2*π*f_amp*t[i]) + 1)
            sig_amp.append( aux1 * aux2 ) ##??????
       
        return np.array(sig_amp)
    
    def noise_signal():
        '''Create noise signals'''
        return rn.uniform(0, 1, len(sig_ph))

    def pink_noise_signal():
        '''Create pink noise signals'''
        return gcn.pink(len(sig_ph), state=None)

    def directionality_signal(delay = 100, white = True, pink = True):
        '''
        Impose directionality, delay one signal X units respect to the other
        
        :param *kwargs:
            delay: delay in system units
            white: sum white noise to the signal
            pink: sum pink noise to the signal
        '''

        end = len(sig_ph) - 1 
        no_coupling_sig = sig_ph + sig_amp

        sig_ph_del = sig_ph[delay:end]
        sig_amp_del = sig_amp[:end-delay+1]
        amp_to_phase_sig = sig_ph_del + sig_amp_del

        sig_ph_del2 = sig_ph[:end-delay+1]
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

        return no_coupling_sig , amp_to_phase_sig, phas_to_amp_sig
        
    sig_ph = phase_signal()
    sig_amp = amplitude_signal()
    sig_white_n = noise_signal()
    sig_pink_n = pink_noise_signal()
    no_coupling_sig , amp_to_phase_sig, phas_to_amp_sig = directionality_signal()


'''compute sinusoidal signal'''
#frequencies = [ 1/20., 1/100, 1/6] #freq, they shall be below one
#amplitudes = [0.5, 1, 2]
#t = np.arange(600) 
#sampling_dt = 1

#sig = create_signal(frequencies, amplitudes, t, gauss = False )#

'''compute delay signal'''
#num_segments = 3
#base_length = 200

#create_victor_signal(num_segments,  base_length)