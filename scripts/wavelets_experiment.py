from numpy.core.numeric import base_repr
from numpy.testing._private.utils import nulp_diff
import pandas as pd
import matplotlib.pyplot as plt
import read_surrogate_data as rsd
import generate_colored_noises as gcn
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



def create_signal(frequencies, amplitudes, t, noise = False, gauss = False):
    ''' return signals as the sum of cosines of different frequencies
    
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

class wavelets_scaling:

    def __init__(self, base, frequency_spacing, num_bands, fourier_factor):
        self.base = base
        self.frequency_spacing = frequency_spacing
        self.num_bands = num_bands
        self.fourier_factor = fourier_factor
   

    
def scales_fourier(wav_kernel, base = 3, frequency_spacing = 0.001, num_bands = 4, fourier_factor = 0.25 ):
    '''
    provides automatic scaling for wavelet spectral frequencies in log2 spacing 
    by using the last fourier mode as initial seed
    :param wav_kernel: provide the wavelet kernel in order to set the scale
    :param *kwargs:
        base: the base over which scan the frequency bands, the bigger, the more spread the bands are
        frequency_spacing: the added spacing in log base scale between frequency bands
        num_bands: number of scales minus one do the spectral decomposition; 
                   ranges form s0  to  s0 * 2^(num_bands+dj)
        fourier_factor: a multiplicative factor for the seed mode, recommended to be 1/4 
    '''
    
    scales = []
    aux = 1/freq_specrum[  np.where(fft1d/len(t) > 0.05)[0][-1] ] #last frequency of  the highest fourier mode    
    s0 = (fourier_factor * aux * ( 1.0 / sampling_dt)) / wav_kernel # 
    #int(np.where(fft1d == max(fft1d[0:nyquist]))[0]/2 )/len(t))
    for i in range(num_bands):
        scales.append(  s0 * base**((num_bands-i)+frequency_spacing) )
    
    return scales

def pywt_compute_wavelets(freq = True):
    
    wavelet = 'cmor1.5-1.0' #kind of wavelet kernel
    if freq:
        scales = (
            np.array(1./np.array(frequencies))
            * (1.0 / sampling_dt)
            * pywt.central_frequency(wavelet)
        )
    else: scales = scales_fourier(pywt.central_frequency(wavelet))
    

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
    return coeffs, 1./np.array(scales)



def niko_compute_wavelets(freq = True):

    k0 = 6 #defines the size of the wavelet kernel, the bigger the smother, but eats up the edges of the data
    wavelet = wa.MorletWavelet()
    central_periods = 1./np.array(frequencies)
    
    if freq: scales = (central_periods * (1.0 / sampling_dt)) / wavelet.fourier_factor(k0)#(
    else: scales = scales_fourier(wavelet.fourier_factor(k0))
        
    
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
    
    return np.array(waves), periods, 1./np.array(scales), cois
    

def amplitude_phase_wav(waves):
    
    amp = []
    phase = []
    for w in waves:
        amp.append( np.abs(waves[0]) )#np.sqrt( coef[0].imag**2 + coef[0].real**2 )
        phase.append(np.angle(waves[0]) )#np.arctan2( coef[0].imag, coef[0].real ) 

    return amp, phase

def wav_reconstructed_signal(sig, waves, no_amp = True, individual_amp = False, global_amp = False ):
    
    '''
    reconstructs the signal in three possible ways, no reescale, rescaling the whole signal, recale individual amplitudes
    
    :param waves: the wavelets that the signal has been decomposed into
    :param **kargs: 
        :no_amplitude: DEFAULT do not reescale the amplitude of the signal, 
                    if one of the other kargs is made True this shall be false
        :individual_amp: if true, reescale each amplitude of each wavelet and then multiply them to reconstruct
        :global_amp: reescale teh whole reconstructed signal to have the same amplitude as the original
    :return: reconstructed signal
    '''
    
    rec_signal = np.zeros(len(waves[0]))  
    slopes_vector = np.ones(len(waves))

    def reconstruction(rec_signal, slopes_vector):
        #recosntruc signal
        for w, m in zip(waves, slopes_vector ):
            rec_signal += np.real(m)*np.abs(w) * np.cos( np.angle(w) ) 
        return rec_signal
        
    def lin_reg_rescaling(rec_signal, slopes_vector):
        rec_signal = reconstruction(rec_signal, slopes_vector)
        #rescaling factor with linear regression
        fit_x = np.vstack([rec_signal, np.ones((rec_signal.shape[0]))]).T
        m, c = np.linalg.lstsq(fit_x, sig)[0]

        return m*rec_signal
    
    def individual_amp_rescaling(rec_signal, slopes_vector):
        #reescale each amplitude and recosntruc signal
        for i, w in enumerate(waves):
            fit_x = np.vstack( [w, np.ones((w.shape[0]))] ).T
            slopes_vector[i], c = np.linalg.lstsq(fit_x, sig)[0]
        rec_signal = reconstruction(rec_signal, slopes_vector) 
        return rec_signal

    if no_amp == individual_amp and no_amp == True:
        raise Exception('one too much True value, provably you forgot to make no_amp False')
    elif no_amp == global_amp and no_amp == True:
        raise Exception('one too much True value, provably you forgot to make no_amp False')
    elif individual_amp == global_amp and global_amp == True:
        raise Exception('what do you pretend making two input argument true!? make one false please')

    if no_amp: return reconstruction(rec_signal, slopes_vector)
    if individual_amp: return individual_amp_rescaling(rec_signal, slopes_vector)
    if global_amp: return lin_reg_rescaling(rec_signal, slopes_vector)


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

    # plotting the fft of the signal in frequency domain
    ax[2].set_title("fft Spectrum of the Signal")
    #ax[2].loglog(fft1d[0:nyquist], 'r')
    ax[2].plot(freq_specrum, fft1d/len(t) )
    ax[2].scatter(frequencies, np.ones( len(frequencies) ) )
    
    # plotting the fft of the signal in scale/length domain and log scale
    ax2 = ax[2].twiny()
    ax2.plot(1./freq_specrum, fft1d/len(t), c='r')
    ax2.scatter(1./np.array(freq_bands_niko), np.ones( len(freq_bands_niko) ), c='r') 
    #ax2.semilogx(len(t)/freq_specrum**2, fft1d/len(t), c='r')
    ax2.set_xscale('log')
    ax2.set_xlabel('wavenum', color='r')
    #ax[2].loglog(frequencies, 'b')
    plt.tight_layout()

def plot_waves_amplitude_phase_WL(title, sig, rec_sig, waves, frequencies):
    
    fig, axs =  plt.subplots(nrows=3, ncols=len(waves),  sharey="row", figsize=(15, 10))#sharex=True,
    fig.suptitle(title, fontsize=20)

    axs[0, 0].set_ylabel("signal")
    axs[1, 0].set_ylabel("amplitude")
    axs[2, 0].set_ylabel("phase")
    
    for i in range(len(waves)):
        axs[0, i].set_title(' F  {0:.2f}  P {1:.2f}'.format(frequencies[i], 1./frequencies[i]) )
        axs[0, i].plot(sig[0: int( len(waves[i])/(i+1) ) ], label = 'original')
        axs[0, i].plot(rec_sig[0: int( len(waves[i])/(i+1) ) ], label = 'reconstructed')
        axs[1, i].plot(np.real(waves[i, :]))
        axs[1, i].plot(np.abs(waves[i, :]))
        axs[2, i].plot(np.angle(waves[i, :]))
    axs[0, 0].legend(loc='upper right', fontsize='small', frameon = False)
    

def plot_comparison_methods(wav1, wav2):
    
    fig2, ax = plt.subplots(4,1)
    # plotting the amplitude
    ax[0].set_xlim(20,len(wav1)-20)
    ymax = max(wav1[0][20:-20]) + 0.3*max(wav2[0][20:-20])
    #ax[0].set_ylim( -ymax, ymax   )
    ax[0].set_title("real wavelets")
    ax[0].plot(wav1.real)
    #ax[0].plot(wav1[0].imag)
    #ax[0].plot(wav2[0].imag)
    ax[0].plot(wav2.real) 

    # plotting the phase
    ax[1].set_xlim(20,len(wav1)-20)
    ax[1].set_title(" phase whavelets")
 
    ax[1].plot(np.angle(wav1))
    ax[1].plot(np.angle(wav2))

    # plotting the amp
    ax[2].set_xlim(20,len(wav1)-20)
    ax[2].set_title(" amp whavelets")
    ax[2].plot(np.abs(wav1))
    ax[2].plot(np.abs(wav2))

    # plotting the reconstruction
    ax[3].set_xlim(20,180)
    ax[3].plot(sig, color ='green')
    ax[3].plot(rec_signal1, color ='blue')
    ax[3].plot(rec_signal2, color ='red')


def FFT(t, sig):
    '''
    compute the fast fourier transform and the frequency scale of the transform
    
    :param t: ordered temporal or spatial sequence
    :param sig: kind of signal to be decomposed
    :return:
        :param f: frequency scale
        :param Y:  fourier tranform of the signal
    '''
    n = len(t)
    Δ = (max(t) - min(t)) / (n-1)
    k = int(n/2)
    freq_escale = np.arange(k) / (n*Δ)
    Y = abs(fft(sig))[:k]
    return (freq_escale, Y)



'''compute signal'''
frequencies = [ 1/20., 1/100, 1/6] #freq, they shall be below one
amplitudes = [0.5, 1, 2]
t = np.arange(600) 
sampling_dt = 1

sig = create_signal(frequencies, amplitudes, t, gauss = True )#


'''compute 1d fourier transformation'''
#fft, ifft
freq_specrum, fft1d = FFT(t, sig)#fft(sig)/len(t)
nyquist = int(len(fft1d))

'''compute wavelet decomposition for 2 different methods'''
waves_pywt, freq_bands_pywt = pywt_compute_wavelets( freq = False )
waves_niko, periods, freq_bands_niko, cois = niko_compute_wavelets(freq = False)

'''reconstruct the signal form the wavelets'''
rec_signal_pywt = wav_reconstructed_signal(sig, waves_pywt, no_amp=False, individual_amp=True)
rec_signal_niko = wav_reconstructed_signal(sig, waves_niko, no_amp=False, individual_amp=True)
#amp, phase = amplitude_phase_wav(coeffs)
#amp2, phase2 = amplitude_phase_wav(waves)

'''plot signals and wavelets'''
plot_signal_phase_fft()
plot_waves_amplitude_phase_WL('python wavelet', sig, rec_signal_pywt, waves_pywt, freq_bands_pywt )
plot_waves_amplitude_phase_WL('niko wavelet', sig, rec_signal_niko, waves_niko, freq_bands_niko)
#plot_comparison_methods()
plt.show()
