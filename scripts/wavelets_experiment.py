import pandas as pd
import matplotlib.pyplot as plt
import read_surrogate_data as rsd
import wavelet_analysis as wa
import numpy as np
from numpy.fft import fft, ifft
import pywt
from scipy import signal

name = '../../package_CMI_prague/data/exp_raw/binfiles/Rossler_bin_0.000.bin'


#df = pd.read_csv('EURUSD.csv',sep='\t', index_col='Date')
#df = rsd.read_bin_bin_dataframe(name)
#df.sort_index(inplace=True)
#df = df.resample('W').last()
#sig =  np.array(df['x'][0:1000])
frequencies = [1/10., 1/100] #items shall be below one
amplitudes = [0.5, 1]
t = np.arange(660) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dt = 1
##dt = 1/len(t)#len(sig)



def create_signal(frequencies, amplitudes, t, noise = False, exp = False):
    ''' return signals as the sume of cosinus of different frequencies
    
    noise or exponential delay can be added default is no'''
    
    sig_cos = []
    sig_exp =  np.real(np.exp(-7*(t-0.4)**2)*np.exp(1j*2*np.pi*2*(t-0.4)))
    sig_noise = np.random.normal(0,0.1, size=(len(t)))
    sig = np.zeros(len(t))

    for i, f in enumerate(frequencies):
        sig_cos.append( np.cos( 2*np.pi * f * t) )  #+ signal.gausspulse(t - 0.4, fc=2)
        sig += amplitudes[i]*sig_cos[i]
    
    if noise: sig += sig_noise
    if exp: sig += sig_exp
    
    return sig
    
sig = create_signal(frequencies, amplitudes, t, noise = True)


#fft, ifft

fft1d = fft(sig)

nyquist = int(len(fft1d)/2.)

def plot_signal_phase_fft():
    plt.rcParams['figure.figsize'] = [12, 7]
    fig, ax = plt.subplots(3,1)
    
    # plotting the signal 
    ax[0].set_xlabel('Time')
    ax[0].set_ylabel('Amplitude')
    ax[0].set_title("Signal")

    ax[0].plot(sig, color ='green')
    
    # plotting the phase spectrum of the signal 
    ax[1].set_title("Phase Spectrum of the Signal")
    ax[1].phase_spectrum(sig, color ='green')

    # plotting the fft of the signal 
    ax[2].set_title("fft Spectrum of the Signal")
    ax[2].plot(fft1d[0:nyquist], 'r')

   
def compute_wavelet():
    #central_periods = int( np.where(fft1d == max(fft1d[0:nyquist]))[0]/2 )
    units = 1
    wavelet = 'cmor1.5-1.0'
    scales = (
                np.array(1./np.array(frequencies))
                * (1.0 / sampling_dt)
                * pywt.central_frequency(wavelet)
            )
    #scales = int(central_periods)

    padded_data =  sig
    coeffs, _ = pywt.cwt(
                padded_data,
                scales=scales,
                wavelet=wavelet,
                sampling_period=1,
                axis=0,
            )
    return coeffs


def amplitude_phase_wab(coef):
    rec_signal = np.zeros(len(coef[0]))
    amp = []
    phase = []
    
    for i, c in enumerate(coef):
        amp.append(np.abs(c))     #np.sqrt( coef[0].imag**2 + coef[0].real**2 )
        phase.append(np.angle(c)) #np.arctan2( coef[0].imag, coef[0].real ) 
        rec_signal += amp[i]*np.cos(phase[i]) #!!!!!!!!!!!!!!!!!!!!!!!!!"""
    
    return amp, phase, rec_signal

sampling_dt = 1.0
k0 = 6
wavelet = wa.MorletWavelet()
central_period = 1./(np.array(frequencies))
scales = (central_period    * (1.0 / sampling_dt)) / wavelet.fourier_factor(k0)#(
           # np.array(1.0/frequency)
          #  * (1.0 / sampling_dt)/wavelet.fourier_factor(k0)
       # )

wave, period, scale, coi = wa.continous_wavelet( 
    sig, dt, pad=True, wavelet=wa.MorletWavelet(), dj =0,  s0 =scales[0] , j1 =0, k0=k0
    )
coeffs = compute_wavelet()
print (coeffs.shape, 'shape coeff')
amp, phase, rec_signal = amplitude_phase_wab(coeffs)
amp2, phase2, rec_signal2 = amplitude_phase_wab(wave)

#print(wave[0])

def plot_amplitude_phase_WL():
    
    fig2, ax = plt.subplots(4,1)
    # plotting the amplitude
    ax[0].set_xlim(20,180)
    ymax = max(coeffs[0][20:-20]) + 0.3*max(coeffs[0][20:-20])
    #ax[0].set_ylim( -ymax, ymax   )
    ax[0].set_title("real im wavelets")
    ax[0].plot(coeffs[0].real)
    ax[0].plot(coeffs[0].imag)
    ax[0].plot(wave[0].imag)
    ax[0].plot(wave[0].real) #!!!!!!!!!!!!!!!!!!!!!!!!!"""

    # plotting the phase
    ax[1].set_xlim(20,180)
    ax[1].set_title(" phase whavelets")
 
    ax[1].plot(phase)
    ax[1].plot(phase2)

    # plotting the amp
    ax[2].set_xlim(20,180)
    ax[2].set_title(" amp whavelets")
    ax[2].plot(amp)
    ax[2].plot(amp2)

    # plotting the reconstruction
    ax[3].set_xlim(20,180)
    ax[3].plot(sig, color ='green')
    ax[3].plot(rec_signal, color ='blue')
    ax[3].plot(rec_signal2, color ='red')


def plot_one_amplitude_phase_WL(coef,  amp, phase, rec_signal2):
    
    fig2, ax = plt.subplots(4,1)
    # plotting the amplitude
    ymax = max(coeffs[0][20:-20]) + 0.3*max(coeffs[0][20:-20])
    #ax[0].set_ylim( -ymax, ymax   )
    ax[0].set_title("real im wavelets")
    ax[0].plot(coef.real)
    ax[0].plot(coef.imag)

    # plotting the phase
    ax[1].set_xlim(20,180)
    ax[1].set_title(" phase whavelets")
 
    ax[1].plot(phase)

    # plotting the amp
    ax[2].set_title("amp whavelets")
    ax[2].plot(amp)

    # plotting the reconstruction
    ax[3].set_xlim(20,180)
    ax[3].plot(sig, color ='green')
    ax[3].plot(rec_signal, color ='blue')

def plot_sig_amp_phase(sig, wave):
    _, axs =  plt.subplots(nrows=4, ncols=2, sharex=True, sharey="row", figsize=(15, 10))
    
    axs[0, 0].set_ylabel("signal")
    axs[1, 0].set_ylabel("Re[wave]")
    axs[2, 0].set_ylabel("phase")
    axs[3, 0].set_ylabel("amplitude")
    
    print ('wave, shape', wave.shape)
    for i in range(wave.shape[0]):
        axs[0, i].set_title(period[i])
        axs[0, i].plot(sig)
        axs[1, i].plot(np.real(wave[i, :]))
        axs[2, i].plot(np.angle(wave[i, :]))
        axs[3, i].plot(np.abs(wave[i, :]))



#print(coeffs)
plot_signal_phase_fft()
#plot_sig_amp_phase(sig, wave)
plot_one_amplitude_phase_WL(coeffs[0], amp[0], phase[0], rec_signal)
plot_one_amplitude_phase_WL(coeffs[1], amp[1], phase[1], rec_signal)
plot_one_amplitude_phase_WL(wave[0], np.abs(wave[0, :]), np.angle(wave[0, :]), rec_signal2)
plt.show()
