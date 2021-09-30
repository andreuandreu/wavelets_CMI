import pandas as pd
import matplotlib.pyplot as plt
import read_surrogate_data as rsd
import wavelet_analysis as wa
import numpy as np
from numpy.fft import fft, ifft
import pywt
from scipy import signal

name = '../data/exp_raw/binfiles/Rossler_bin_0.000.bin'


#df = pd.read_csv('EURUSD.csv',sep='\t', index_col='Date')
#df = rsd.read_bin_bin_dataframe(name)
#df.sort_index(inplace=True)
#df = df.resample('W').last()
#sig =  np.array(df['x'][0:1000])
frequency = 1/10.
t = np.arange(200) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##dt = 1/len(t)#len(sig)
#sig  = np.cos(2 * np.pi * 7 * t) + np.real(np.exp(-7*(t-0.4)**2)*np.exp(1j*2*np.pi*2*(t-0.4)))
sig  = np.cos( 2*np.pi * frequency * t) + np.random.normal(0,0.1, size=(len(t))) #+ signal.gausspulse(t - 0.4, fc=2)
#sig  = np.cos(2*np.pi*1*t) + np.cos(2*np.pi*10*t) + np.cos(2*np.pi*100*t)



#fft, ifft

fft1d = fft(sig)

nyquist=int(len(fft1d)/2.)

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

    #ploting the wavelwet transform
    #ax[2].pcolormesh(t, freq, np.abs(cwtmatr), cmap='viridis', shading='gouraud')
    #ax[2].imshow(cwtmatr, extent=[-1, 1, 1, 31], cmap='PRGn', aspect='auto',
    #            vmax=abs(cwtmatr).max(), vmin=-abs(cwtmatr).max())
 


#phase  = np.arctan()


dt = 1
def compute_wavelet():
    central_periods = int( np.where(fft1d == max(fft1d[0:nyquist]))[0]/2 )
    print('maxxxxxx', central_periods )
    units = 1
    wavelet = 'cmor1.5-1.0'
    scales = (
                np.array(central_periods)
                * (1.0 / dt)
                * pywt.central_frequency(wavelet)
            )
    scales = int(central_periods)
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
    
    amp = np.abs(coef[0])#np.sqrt( coef[0].imag**2 + coef[0].real**2 )
    print( 'amp', amp[0])
    phase = np.angle(coef[0])#np.arctan2( coef[0].imag, coef[0].real ) 
    rec_signal = amp*np.cos(phase) #!!!!!!!!!!!!!!!!!!!!!!!!!"""
    return amp, phase, rec_signal

sampling_dt = 1.0
k0 = 6
wavelet = wa.MorletWavelet()
central_period = 1./frequency
scales = (central_period    * (1.0 / sampling_dt)) / wavelet.fourier_factor(k0)#(
           # np.array(1.0/frequency)
          #  * (1.0 / sampling_dt)/wavelet.fourier_factor(k0)
       # )

wave, period, scale, coi = wa.continous_wavelet( 
    sig, dt, pad=True, wavelet=wa.MorletWavelet(), dj =0,  s0 =scales , j1 =0, k0=k0
    )
coeffs = compute_wavelet()

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



#print(coeffs)
plot_signal_phase_fft()
plot_amplitude_phase_WL()
plt.show()
