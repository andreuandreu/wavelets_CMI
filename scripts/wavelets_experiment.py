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
frequencies = [1/20., 1/100, 1/6] #items shall be below one
amplitudes = [0.5, 1, 2]
t = np.arange(400) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
sampling_dt = 1
frequency_spacing = 1 #the spacing in log 2 scale
num_bands = 11 # number of scales minus one do the spectral decomposition; ranges form s0  to  s0 * 2^(num_bands+dj)` 

##dt = 1/len(t)#len(sig)


def create_signal(frequencies, amplitudes, t, noise = False, gauss = False):
    ''' return signals as the sume of cosinus of different frequencies
    
    noise or exponential delay can be added default is no'''
    
    sig_cos = []
    sig_gauss =  10*np.real(np.exp(-0.001*(t-len(t)/2)**2))#*np.exp(1j*2*np.pi*2*(t-0.4)))
    sig_noise = np.random.normal(0,0.1, size=(len(t)))
    sig = np.zeros(len(t))

    for i, f in enumerate(frequencies):
        sig_cos.append( np.cos( 2*np.pi * f * t) )  #+ signal.gausspulse(t - 0.4, fc=2)
        sig += amplitudes[i]*sig_cos[i]
    
    if noise: sig += sig_noise
    if gauss: sig += sig_gauss
    
    return sig
    



 





def pywt_compute_wavelets():
    

    units = 1
    wavelet = 'cmor1.5-1.0'
    scales = (
                np.array(1./np.array(frequencies))
                * (1.0 / sampling_dt)
                * pywt.central_frequency(wavelet)
            )
    #scales = int(central_periods)
    #scales = central_periods
    padded_data =  sig
    coeffs, _ = pywt.cwt(
                padded_data,
                scales=scales,
                wavelet=wavelet,
                sampling_period=1,
                #axis=0,
            )
    return coeffs



def niko_compute_wavelets(freq = True):

    k0 = 6 #defines the size of the wavelet kernel, the bigger the smother, but eats up the edges of the data
    wavelet = wa.MorletWavelet()
    central_periods = 1./np.array(frequencies)
    
    if freq: scales = (central_periods    * (1.0 / sampling_dt)) / wavelet.fourier_factor(k0)#(
    else: 
        scales = []
        s0 = int( np.where(fft1d == max(fft1d[0:nyquist]))[0]/2 )/len(t)
        for i in enumerate(num_bands):
            scales.append(  s0 * 2**(i+frequency_spacing) )
    
    waves = []
    periods = []
    wav_scales = []
    cois = []
    
    for s in scales:
        wave, period, scale, coi = wa.continous_wavelet( 
            sig, sampling_dt, pad=True, wavelet=wa.MorletWavelet(), dj =0,  s0 =s , j1 =0, k0=k0
            )
        waves.append(wave[0])
        periods.append(period)
        wav_scales.append(scale)
        cois.append(coi)

    return np.array(waves), periods, scales, cois
    

def amplitude_phase_wav(waves):
    
    amp = []
    phase = []
    for w in waves:
        amp.append( np.abs(waves[0]) )#np.sqrt( coef[0].imag**2 + coef[0].real**2 )
        phase.append(np.angle(waves[0]) )#np.arctan2( coef[0].imag, coef[0].real ) 

    return amp, phase

def wav_reconstructed_signal(waves):
    
    rec_signal = np.zeros(len(waves[0]))

    #recosntruc signal  
    for w in waves:
        rec_signal += np.abs(w) * np.cos( np.angle(w) ) 
   
    #rescaling factor with linear regression
    fit_x = np.vstack([rec_signal, np.ones((rec_signal.shape[0]))]).T
    m, c = np.linalg.lstsq(fit_x, sig)[0]
    
    return m*rec_signal


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
    ax[2].loglog(fft1d[0:nyquist], 'r')
    ax[2].loglog(frequencies, 'b')



def plot_waves_amplitude_phase_WL(title, sig, rec_sig, waves):
    
    fig, axs =  plt.subplots(nrows=3, ncols=len(waves), sharex=True, sharey="row", figsize=(15, 10))
    fig.suptitle(title, fontsize=20)

    axs[0, 0].set_ylabel("signal")
    axs[1, 0].set_ylabel("ampitude")
    axs[2, 0].set_ylabel("phase")
    
    for i in range(len(waves)):
        axs[0, i].set_title(frequencies[i])
        axs[0, i].plot(sig, label = 'original')
        axs[0, i].plot(rec_sig, label = 'reconstructed')
        axs[1, i].plot(np.real(waves[i, :]))
        axs[1, i].plot(np.abs(waves[i, :]))
        axs[2, i].plot(np.angle(waves[i, :]))
    axs[0, 0].legend(loc='upper right', fontsize='small', frameon = False)
 

def plot_comparison_methods(wav1, wav2):
    
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


def FFT(t, sig):
    n = len(t)
    Δ = (max(t) - min(t)) / (n-1)
    k = int(n/2)
    f = np.arange(k) / (n*Δ)
    Y = abs(fft(sig))[:k]
    return (f, Y)

#compute signal
sig = create_signal(frequencies, amplitudes, t )#gauss = True

#compute 1d fourier trasformation
#fft, ifft
freq_specrum, fft1d = FFT(t, sig)#fft(sig)/len(t)
nyquist = int(len(fft1d)/2.)

#compute wavelet decomposition for 2 different methods
waves_pywt = pywt_compute_wavelets()
waves_nico, periods, scales, cois = niko_compute_wavelets()


rec_signal_pywt = wav_reconstructed_signal(waves_pywt)
rec_signal_nico = wav_reconstructed_signal(waves_nico)
#amp, phase = amplitude_phase_wav(coeffs)
#amp2, phase2 = amplitude_phase_wav(waves)

#plot signals and wavelets
plot_signal_phase_fft()
#plot_waves_amplitude_phase_WL('python wavelet', sig, rec_signal_pywt, waves_pywt)
#plot_waves_amplitude_phase_WL('nico wavelet', sig, rec_signal_nico, waves_nico)
#plot_amplitude_phase_WL()
plt.show()
