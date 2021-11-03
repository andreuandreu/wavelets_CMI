from numpy.core.numeric import base_repr
from numpy.testing._private.utils import nulp_diff
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import read_surrogate_data as rsd

import wavelet_analysis as wa
import numpy as np
from numpy.fft import fft, ifft
import pywt
from scipy import signal
import create_synthetic_signals as css


'''
Functions to decompose a signal into its component continuos wavelets and reconstruct it


- 1st provide a univaluated signal, chose one function from the create_synthetic_signals script

- 2nd do the fast fourier transfrom of the signal using the FFT(t, sig) function

- 3rd compute wavelet decomposition for one or both provided methods
    niko_compute_wavelets(freq = True)
    pywt_compute_wavelets(freq = True)
    if frequency is true, then the decomposition would be around the provided list of frequencies
    otherwise it will compute a discrete set of frequencies in base 2 to decompose
    in this case one has to initialize the values of the desired intervals into the wavelets_scaling class

- 4th reconstruct the signal using function  wav_reconstructed_signal(sig, waves_pywt, no_amp=False, individual_amp=True)
    here 3 methods of reconstruction can be probided that rescale the amplitudes, 
        no_amplitude: does not rescale
        global_amp: rescales the sum of the wavelets to the amplitude of the signal
        individual_amp: rescales each wavelet to the amplitude of the signal

- 5th plots
        plot_signal_phase_fft(time, signal): signal + phase transform + FFt in frequency and wavenumber domain
        plot_waves_amplitude_phase_WL(): reconstruction of the signal and each wavelet (phase and amplitude) corresponding to the provided frequency
        plot_comparison_methods(): to compare the wavelets and reconstructions of two different wavelet methods ->in development<-   
'''


class wavelets_scaling:
    '''
    Initialize the parameters to be fed in the "automatic" band selection for the wavelet decomposition in 
    regular bands.

    :param base: the base over which scan the frequency bands, the bigger, the more spread the bands are
    :param frequency_spacing: the added spacing in log base scale between frequency bands
    :param num_bands: number of scales minus one do the spectral decomposition; 
                   ranges form s0  to  s0 * 2^(num_bands+dj)
    :param fourier_factor: a multiplicative factor for the seed mode, recommended to be 1/4 
    '''
    def __init__(self, base = 2, frequency_spacing = 0.001, num_bands = 4, fourier_factor = 0.25):
        self.base = base
        self.frequency_spacing = frequency_spacing
        self.num_bands = num_bands
        self.fourier_factor = fourier_factor
   

    
def scales_fourier(wav_kernel, par ): #base = 3, frequengit branch -dcy_spacing = 0.001, num_bands = 4, fourier_factor = 0.25
    '''
    provides automatic scaling for wavelet spectral frequencies in log2 spacing 
    by using the last fourier mode as initial seed
    :param wav_kernel: provide the wavelet kernel in order to set the scale
    :param par: class containing the parameters to be used (base, frequency_spacing, num_bands_ fourier_factor)
    :type par: class with self initializing parameter values
    '''
    
    print(par.frequency_spacing, par.base, par.num_bands)
    scales = []
    aux = 1/freq_spectrum[  np.where(fft1d/len(t) > 0.05)[0][-1] ] #last frequency of  the highest fourier mode    
    s0 = (par.fourier_factor * aux * ( 1.0 / sampling_dt)) / wav_kernel # 
    #int(np.where(fft1d == max(fft1d[0:nyquist]))[0]/2 )/len(t))
    for i in range(par.num_bands):
        scales.append(  s0 * par.base**((par.num_bands-i)+par.frequency_spacing) )
    
    return scales

def pywt_compute_wavelets(kernel_pywl = 'cmor', freq = True, par_scales =  wavelets_scaling ):
    
    
    if freq:
        scales = 1/frequencies
    else: 
        scales = scales_fourier( pywt.central_frequency(kernel_pywl), par_scales)
    

    #scales = int(central_periods)
    
    padded_data =  sig
    coeffs, freq_pywt = pywt.cwt(
                padded_data,
                scales=scales,
                wavelet=kernel_pywl,
                sampling_period=sampling_dt,
                #axis=0,
            )
    return coeffs, freq_pywt



def niko_compute_wavelets(freq = True, par_scales = wavelets_scaling):

    k0 = 9 #defines the size of the wavelet kernel, the bigger the smother, but eats up the edges of the data
    wavelet = wa.MorletWavelet()
    central_periods = sampling_dt/np.array(frequencies)
    
    if freq: scales = (central_periods ) / wavelet.fourier_factor(k0)#(
    else: scales = scales_fourier(wavelet.fourier_factor(k0), par_scales)
        
    
    waves = []
    wav_periods = []
    wav_scales = []
    cois = []
    
    for s in scales:
        wave, period, scale, coi = wa.continous_wavelet( 
            sig, sampling_dt, pad=True, wavelet=wa.MorletWavelet(), dj =0,  s0 =s , j1 =0, k0=k0
            )
        waves.append(wave[0])
        wav_periods.append(period[0])
        wav_scales.append(scale[0])
        cois.append(coi)
    print ('\n\n this is that', wav_periods,  '\n',wav_scales,'\n' )
    return np.array(waves), np.array(wav_periods), 1./np.array(wav_scales), cois
    

def amplitude_phase_wav(waves, name_file):
    
    amp = []
    phase = []
    for w in waves:
        amp.append( np.abs(w) )#np.sqrt( coef[0].imag**2 + coef[0].real**2 )
        phase.append(np.angle(w) )#np.arctan2( coef[0].imag, coef[0].real ) 

    np.save(name_file+'_amp.npy', amp, allow_pickle=True, fix_imports=True)
    np.save(name_file+'_pha.npy', phase, allow_pickle=True, fix_imports=True)
    return amp, phase

def wav_reconstructed_signal(sig, waves, no_amp = True, individual_amp = False, global_amp = False ):
    
    '''
    reconstructs the signal in three possible ways, no reescale, rescaling the whole signal, recale individual amplitudes
    
    :param waves: the wavelets that the signal has been decomposed into
    :param **kwargs: 
        :no_amplitude: DEFAULT do not rescale the amplitude of the signal, 
                    if one of the other kwargs is made True this shall be false
        :individual_amp: if true, rescale each amplitude of each wavelet and then multiply them to reconstruct
        :global_amp: rescale teh whole reconstructed signal to have the same amplitude as the original
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


def plot_signal_plus_average(ax, time, signal, time_ave, signal_ave, average_over = 5):

    # plotting the signal 
    ax.set_xlabel('Time')
    ax.set_ylabel('Amplitude')
    ax.set_title("Signal + Time Average")
    ax.set_xlim([t[0], t[-1]])
    ax.plot(t, sig, color ='green')

    #plotting the average of the signal
    ax.plot(time_ave, signal_ave, label = 'time average (n={})'.format(6))
    ax.legend(loc='upper right', fontsize='small', frameon = False)



def plot_signal_phase_fft(time, signal, freq_spectrum, frequencies, fft1d):

    plt.rcParams['figure.figsize'] = [12, 7]
    fig, ax = plt.subplots(3,1)
    
    # plotting the signal 
    ax[0].set_xlabel('Time')
    ax[0].set_ylabel('Amplitude')
    ax[0].set_title("Signal + Time Average")
    ax[0].set_xlim([time[0], time[-1]])
    ax[0].plot(time, signal, color ='green')
    ax[0].plot(t, sig, color ='r')
    #plotting the average of the signal
    time_ave, signal_ave = css.get_ave_values(time, signal, 6)
    #ax[0].plot(time_ave, signal_ave, label = 'time average (n={})'.format(6))
    ax[0].legend(loc='upper right', fontsize='small', frameon = False)
   
    # plotting the phase spectrum of the signal 
    ax[1].set_title("Phase Spectrum of the Signal")
    ax[1].phase_spectrum(signal, color ='green')

    # plotting the fft of the signal in frequency domain
    ax[2].set_title("fft Spectrum of the Signal")
    #ax[2].loglog(fft1d[0:nyquist], 'r')
    ax[2].plot(freq_spectrum, 2*fft1d/len(t) )
    ax[2].scatter(frequencies, np.ones( len(frequencies) ) )
    
    # plotting the fft of the signal in scale/length domain and log scale
    ax2 = ax[2].twiny()
    ax2.plot(1./freq_spectrum, fft1d/len(time), c='r')
    #ax2.scatter(1./np.array(freq_bands_niko), np.ones( len(freq_bands_niko) ), c='r') 
    #ax2.semilogx(len(t)/freq_spectrum**2, fft1d/len(t), c='r')
    ax2.set_xscale('log')
    ax2.set_xlabel('wavenum', color='r')
    #ax[2].loglog(frequencies, 'b')
    plt.tight_layout()



def plot_fft(ax, plot_direction='vertical', yticks=None, ylim=None):

    if plot_direction == 'vertical':


        #fft with scales (periods/dimension) in log scale
        scales = 1./freq_spectrum
        scales_log = np.log2(scales)
        ax.scatter( [np.mean(2*fft1d/len(t))]* len(frequencies) , np.log2(1/frequencies))
        ax.plot( 2*fft1d/len(t), scales_log, 'r-', label='Fourier Transform')
        ax.legend(loc='upper right', fontsize='small', frameon = False)
        #plot power
        #variance = np.std(sig)**2
        #ax.plot( 2* variance * abs(fft1d) ** 2/len(fft1d), scales_log, 'k--', linewidth=1, label='FFT Power Spectrum')
        
        #set log scale
        ax.set_yticks(np.log2(yticks))
        ax.set_yticklabels(yticks)
        ax.invert_yaxis()
        ax.set_ylim(ylim[0], -1)

        #fft with frequencies and linear scale
        ax2 = ax.twinx()
        ax2.plot( 2*fft1d/len(t), freq_spectrum, 'g--' )
        ax2.set_yscale('linear')
        ax2.set_ylabel('Frequency', color='g')


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
    ax[3].plot(rec_signal_pywt, color ='blue')
    ax[3].plot(rec_signal_niko, color ='red')



def plot_scalogram_fft_signal_together(time, signal, time_ave, signal_ave, freq, wavelets, waveletname):
   
    #prepare plot
    fig = plt.figure(figsize=(12,12))
    spec = gridspec.GridSpec(ncols=6, nrows=6)
    top_ax = fig.add_subplot(spec[0, 0:5])
    bottom_left_ax = fig.add_subplot(spec[1:, 0:5])
    bottom_right_ax = fig.add_subplot(spec[1:, 5])

    #plot signal
    plot_signal_plus_average(top_ax, time, signal, time_ave, signal_ave, average_over = 5)
    
    #plot scalogram
    yticks, ylim, im = plot_wavelet_scalogram(time, freq, wavelets, waveletname, bottom_left_ax )
    #position colorbar
    cbar_ax = fig.add_axes([0.9, 0.85, 0.01, 0.1])
    fig.colorbar(im, cax=cbar_ax, orientation="vertical")
    
    #plot fft
    plot_fft(bottom_right_ax, plot_direction='vertical', yticks=yticks, ylim=ylim)
    bottom_right_ax.set_ylabel('Period [years]', fontsize=14)
    plt.tight_layout()
 


def plot_wavelet_scalogram(time, freq, wavelets, waveletname = 'cmor', ax = plt.subplots(figsize=(15, 10))  ):

    #prepare data
    power = (abs(wavelets)) ** 2
    period = 1 / (freq)  # 1/freq
    
    #prepare contours levels 
    levels = [0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8]#[0.0625, 0.25,  1, 8]##[8, 16000, 32000, 128000]#
    contourlevels = np.log2(levels)
    cmap = plt.cm.seismic
    im = ax.contourf(time, np.log2(period), np.log2(power), contourlevels, extend='both',cmap=cmap)
    
    #prepare axis tags
    ax.set_title('Wavelet '+waveletname+' Transform/time Power Spectrum)', fontsize=20)
    ax.set_ylabel('Scale (years)', fontsize=18)
    ax.set_xlabel('Time', fontsize=18)
    
    #arrange y axis scales
    yticks = 2**np.arange(np.ceil(np.log2(period.min())), np.ceil(np.log2(period.max())))
    ax.set_yticks(np.log2(yticks))
    ax.set_yticklabels(yticks)
    ax.invert_yaxis()
    ylim = ax.get_ylim()
    ax.set_ylim(ylim[0], -1)
    
    return yticks, ylim, im



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
seed_freq = [ 1/20., 1/100, 1/6] #freq, they shall be below one
amplitudes = [0.5, 1, 2]

#sampling_dt = 1
#t = np.arange(600) 
#sig = css.create_signal(seed_freq, amplitudes, t, gauss = True, noise = True )#

data = './data/output/'#'../data/' in case you are in the scripts folder
#name_sig = 'ENSO34_1880-2020'
#t, sig = css.online_ENSO_34()


name_ENSO = './data/imput/s_nino3418702020.csv'
name_rain = './data/imput/s_allindiarain18712016.csv'
name_sig = 'rain_india_manuel'
#name_sig = 'ENSO_manuel'
#name_sig = 'ENSO_online'
t, sig = css.read_ENSO_rain_manuel_files(name_rain)
#sig, dates = css.online_ENSO_34()
time_ave, signal_ave = css.get_ave_values(t, sig, 3)
#t, sig = css.get_ave_values(t, sig, 3)

print(sig, t)
sampling_dt = t[1]-t[0]
print(sampling_dt, t[-1], len(t))
#frequencies = 1/np.arange(1, len(sig)//12, 8)
frequencies = 1/np.array([0.083, 0.1,0.5,0.9,1,2,4,5,6,7,8, 16, 32,64, 126])#1/np.arange(1, 128)
#frequencies = np.array([  1., 1.1, 0.25, 0.125, 0.0039 ])/sampling_dt# np.arange(1, 256)*0.0039


'''compute 1d fourier transformation'''
#fft, ifft
freq_spectrum, fft1d = FFT(t, sig)#fft(sig)/len(t)
nyquist = int(len(fft1d))


'''compute wavelet decomposition for 2 different methods'''
kernel_pywl = 'gaussian'#'cmor1.5-1.0'#'cmor'# #kind of wavelet kernel
kernel_niko = 'morlet'
bands_par = wavelets_scaling(num_bands=6)
waves_pywt, freq_bands_pywt = pywt_compute_wavelets( freq = True, par_scales = bands_par )
waves_niko, periods_niko, freq_bands_niko, cois = niko_compute_wavelets(freq = True, par_scales =  bands_par)

fig, ax = plt.subplots(1)
ax.plot(freq_spectrum)
ax.plot(np.arange(0, len(periods_niko))*len(freq_spectrum) /
        len(periods_niko), 1/periods_niko[::-1])
print ('\n periods niko ', periods_niko, '\n')
print('\n periods niko ', 1/periods_niko, '\n')

name_files_pywt = data + name_sig +'_wavelet_vecors_pywt_'+ kernel_pywl
amplitude_phase_wav(waves_pywt, name_files_pywt)

'''reconstruct the signal form the wavelets'''

rec_signal_pywt = wav_reconstructed_signal(sig, waves_pywt, no_amp=False, individual_amp=True)
rec_signal_niko = wav_reconstructed_signal(sig, waves_niko, no_amp=False, individual_amp=True)
#amp, phase = amplitude_phase_wav(coeffs)
#amp2, phase2 = amplitude_phase_wav(waves)

'''plot signals and wavelets'''
plot_signal_phase_fft(t, sig, freq_spectrum, frequencies, fft1d)
plot_scalogram_fft_signal_together(t, sig, time_ave, signal_ave, 1/periods_niko , waves_niko, waveletname = kernel_niko)



#time_ave2, signal_ave2 = css.get_ave_values(t, sig, 3)
t, sig = css.online_ENSO_34()
sampling_dt = t[1]-t[0]
freq_spectrum, fft1d = FFT(t, sig)#fft(sig)/len(t)
waves_niko, periods_niko, freq_bands_niko, cois = niko_compute_wavelets( freq = True, par_scales = bands_par )
plot_signal_phase_fft(time_ave, signal_ave, freq_spectrum, frequencies, fft1d)
plot_scalogram_fft_signal_together(t, sig, time_ave, signal_ave, 1/periods_niko, waves_niko, waveletname = kernel_pywl)



#plot_wavelet_scalogram(t, freq_bands_niko, waves_niko, waveletname = kernel_niko  )
#plot_wavelet_scalogram(t, freq_bands_pywt, waves_pywt, waveletname = kernel_pywl  )
#plot_waves_amplitude_phase_WL('python wavelet', sig, rec_signal_pywt, waves_pywt, freq_bands_pywt )
#plot_waves_amplitude_phase_WL('niko wavelet', sig, rec_signal_niko, waves_niko, freq_bands_niko)
#plot_comparison_methods()
plt.show()