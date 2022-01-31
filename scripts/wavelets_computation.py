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




-  compute wavelet decomposition for one or both provided methods
    pywt_compute_wavelets( )
    niko_compute_wavelets( )
    
-   reconstruct the signal using function  wav_reconstructed_signal(sig, waves_pywt, no_amp=False, individual_amp=True)
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
   

# base = 3, frequengit branch -dcy_spacing = 0.001, num_bands = 4, fourier_factor = 0.25
def scales_fourier(freq_spectrum, fft1d, sampling_dt, wav_kernel, par=wavelets_scaling()):
    '''
    provides automatic scaling for wavelet spectral frequencies in log2 spacing 
    by using the last fourier mode as initial seed
    :param freq_spectrum: the spectrum of frequencies where the wavelet will be computed
    :param fft1d: the Fourier transform of a guiven signal to be analysed 
    :param wav_kernel: provide the wavelet kernel in order to set the scale
    :param par: class containing the parameters to be used (base, frequency_spacing, num_bands_ fourier_factor)
    :type par: class with self initializing parameter values
    '''
    
    print(par.frequency_spacing, par.base, par.num_bands)
    scales = []
    aux = 1/freq_spectrum[  np.where(fft1d/len(fft1d) > 0.05)[0][-1] ] #last frequency of  the highest fourier mode    
    s0 = (par.fourier_factor * aux * ( 1.0 / sampling_dt)) / wav_kernel # 
    #int(np.where(fft1d == max(fft1d[0:nyquist]))[0]/2 )/len(t))
    for i in range(par.num_bands):
        scales.append(  s0 * par.base**((par.num_bands-i)+par.frequency_spacing) )
    
    return scales


def pywt_compute_wavelets(sig, frequencies, kernel_name='cmor1.5-1.0'):
    
    '''
    a warp to compute the wavelet using the pywl package. 

    :param sig: the signal to be analysed
    :param frequencies: the individual frequencies where the wavelet shall be computed,
                        be aware of the units of that!!
    :param kernel_name: the wavelet kernel to be used
    :param freq: 
    '''
    scales = 1/frequencies
    coeffs, freq_pywt = pywt.cwt(
                sig,
                scales=scales,
                wavelet= kernel_name,
                sampling_period=1.0,
                #axis=0,
            )

    #print('peeeriod or', 1/freq_pywt, '\n')
    return coeffs, freq_pywt#*sampling_dt


def niko_compute_wavelets(sig, frequencies, sampling_dt, kernel_name = 'morlet' ):

    k0 = 6. #defines the size of the wavelet kernel, the bigger the smother, but eats up the edges of the data
    if kernel_name == 'morlet':
        wavelet_kernel = wa.MorletWavelet()
    else:
        print('give me a good wavelet kernel name for niko, duh!!')
        exit()
    central_periods = 1./np.array(frequencies)
    
    scales = (central_periods )*(1.0/sampling_dt) / wavelet_kernel.fourier_factor(k0)#(
        
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
    
def read_wavelets(name_base):


    amplitude = np.load(name_base+'_amp.npy', fix_imports=False)
    phase = np.load(name_base+'_pha.npy', fix_imports=False)
    scales = np.load(name_base+'_ska.npy', fix_imports=False)
    return amplitude, phase, scales

def write_amplitude_phase_scale_wav(waves, scales, name_file):
    
    amp = []
    phase = []
    for w in waves:
        amp.append( np.abs(w) )#np.sqrt( coef[0].imag**2 + coef[0].real**2 )
        phase.append(np.angle(w) )#np.arctan2( coef[0].imag, coef[0].real ) 

    np.save(name_file+'_amp.npy', amp, allow_pickle=True, fix_imports=True)
    np.save(name_file+'_pha.npy', phase, allow_pickle=True, fix_imports=True)
    np.save(name_file+'_ska.npy', scales, allow_pickle=True, fix_imports=True)
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
            slopes_vector[i], c = np.linalg.lstsq(fit_x, sig, rcond=None)[0]
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