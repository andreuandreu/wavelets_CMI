from numpy.core.numeric import base_repr
from numpy.testing._private.utils import nulp_diff
import pandas as pd
import matplotlib.pyplot as plt

import wavelet_analysis as wa
import numpy as np
from numpy.fft import fft, ifft
import pywt
from scipy import signal
import provide_signals as css


'''
Functions to decompose a signal into its component continuos wavelets and reconstruct it


- 1st provide a univaluated signal, chose one function from the provide_signals script

- 2nd settle in the units of the signal by chosing a value for sampling_dt

- 3rd do the fast fourier transfrom of the signal using the FFT(t, sig) function

- 4th compute wavelet decomposition for one or both provided methods
    pywt_compute_wavelets( kernel_pywl=kernel_pywl, freq=True,  par_scales=bands_par)
    pywt_compute_wavelets(freq = True)
    if frequency is true, then the decomposition would be around the provided list of frequencies
    otherwise it will compute a discrete set of frequencies in base 2 to decompose
    in this case one has to initialize the values of the desired intervals into the wavelets_scaling class

- 5th reconstruct the signal using function  wav_reconstructed_signal(sig, waves_pywt, no_amp=False, individual_amp=True)
    here 3 methods of reconstruction can be probided that rescale the amplitudes, 
        no_amplitude: does not rescale
        global_amp: rescales the sum of the wavelets to the amplitude of the signal
        individual_amp: rescales each wavelet to the amplitude of the signal


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


def pywt_compute_wavelets(sig, kernel_name='cmor1.5-1.0', freq=True, par_scales=wavelets_scaling):
    
    
    if freq:
        scales = 1/frequencies
    else: 
        scales = scales_fourier( pywt.central_frequency(kernel_name), par_scales)
    
    coeffs, freq_pywt = pywt.cwt(
                sig,
                scales=scales,
                wavelet= kernel_name,
                sampling_period=1.0,
                #axis=0,
            )

    #print('peeeriod or', 1/freq_pywt, '\n')

    return coeffs, freq_pywt#*sampling_dt


def niko_compute_wavelets(sig, kernel_name = 'morlet' , freq = True, par_scales = wavelets_scaling):

    k0 = 9. #defines the size of the wavelet kernel, the bigger the smother, but eats up the edges of the data
    if kernel_name == 'morlet':
        wavelet_kernel = wa.MorletWavelet()
    else:
        print('give me a good wavelet kernel name for niko, duh!!')
        exit()
    central_periods = 1./np.array(frequencies)
    
    if freq: scales = (central_periods )*(1.0/sampling_dt) / wavelet_kernel.fourier_factor(k0)#(
    else: scales = scales_fourier(wavelet_kernel.fourier_factor(k0), par_scales)
        
    
    waves = []
    wav_periods = []
    wav_scales = []
    cois = []
    
    for s in scales:
        wave, period, scale, coi = wa.continous_wavelet( 
            sig, sampling_dt, pad=True, wavelet=wavelet_kernel, dj =0,  s0 =s , j1 =0, k0=k0
            )
        waves.append(wave[0])
        wav_periods.append(period[0]*sampling_dt)
        wav_scales.append(scale[0]*sampling_dt)
        cois.append(coi)

    return np.array(waves), np.array(wav_periods), 1./np.array(wav_scales), cois
    

def write_amplitude_phase_wav(waves, name_file):
    
    amp = []
    phase = []
    for w in waves:
        amp.append( np.abs(w) )#np.sqrt( coef[0].imag**2 + coef[0].real**2 )
        phase.append(np.angle(w) )#np.arctan2( coef[0].imag, coef[0].real ) 

    np.save(name_file+'_amp.npy', amp, allow_pickle=True, fix_imports=True)
    np.save(name_file+'_pha.npy', phase, allow_pickle=True, fix_imports=True)
    print(' \n saving amplitude ', name_file+'_amp.npy')
    print(' saving phase in ', name_file+'_pha.npy\n')
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

'''folder for data output'''
data = './data/output/'  # '../data/' in case you are in the scripts folder

'''compute signal
seed_freq = [1/20., 1/100, 1/6]  # freq, they shall be below one
amplitudes = [0.5, 1, 2]
#sampling_dt = 1
t = np.arange(600)

t, sig = css.create_signal(seed_freq, amplitudes, t, gauss=gauss, noise=noise)
sig_tag = 'synthetic'
'''
gauss = True
noise = True
name_source = {
    'rain_india_manuel': './data/imput/s_allindiarain18712016.csv',
    'ENSO_manuel': './data/imput/s_nino3418702020.csv',
    'ENSO_online': 'http://paos.colorado.edu/research/wavelets/wave_idl/sst_nino3.dat',
    'synthetic': 'gauss_' + str(gauss)[0] + 'noise_' + str(noise)[0]
}

#sig_tag = 'rain_india_manuel'
sig_tag = 'ENSO_manuel'
t, sig = css.read_ENSO_rain_manuel_files(name_source[sig_tag])
#t, sig = css.online_ENSO_34()

if sig_tag == 'ENSO_manuel':
    sig = sig*20

time = np.arange(0, len(sig))
#t, sig = css.get_ave_values(t, sig, 3)

'''characteristics of the signal and the processing'''
unit = 'month'
sampling_dt = 1# t[1]-t[0]
print('\n sampling period, last time and length of time', sampling_dt, time[-1], len(time), '\n')
#frequencies = 1/np.array([0.083, 0.1,0.5,0.9,1,2,4,5,6,7,8, 16, 32,64, 126])#
frequencies = 1/((np.arange(1, 256)*sampling_dt))  #

'''compute 1d fourier transformation'''
#fft, ifft
freq_spectrum, fft1d = FFT(time, sig)#fft(sig)/len(t)
nyquist = int(len(fft1d))


'''compute wavelet decomposition for 2 different methods'''
kernel_pywl = 'cmor1.5-1.0'  # 'cmor'# #kind of wavelet kernel'gaussian'#
kernel_niko = 'morlet'
bands_par = wavelets_scaling(num_bands=6)
waves_pywt, freq_bands_pywt = pywt_compute_wavelets( sig,
    kernel_name=kernel_pywl, freq=True,  par_scales=bands_par)
waves_niko, periods_niko, freq_bands_niko, cois = niko_compute_wavelets( sig,
    kernel_name =  kernel_niko, freq = True, par_scales =  bands_par)


'''satore the amplitude and phase of the waveleets in numpy files'''
name_files_pywt = data + sig_tag +'_'+ unit +'_wavelet_vecors_pywt_'+ kernel_pywl 
write_amplitude_phase_wav(waves_pywt, name_files_pywt)



'''reconstruct the signal form the wavelets'''
rec_signal_pywt = wav_reconstructed_signal(sig, waves_pywt, no_amp=False, individual_amp=True)
rec_signal_niko = wav_reconstructed_signal(sig, waves_niko, no_amp=False, individual_amp=True)