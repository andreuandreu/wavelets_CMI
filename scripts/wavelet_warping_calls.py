from turtle import end_fill
import numpy as np
import provide_signals as ps
import wavelets_computation as wc
import wavelets_ploting as wp
import configparser
import surogates as su
import matplotlib.pyplot as plt
import os
'''
code that warps the functions necessary tu do a preliminary analysis of a guiven signal in 
wavelets, save them in .npy files and plot them.



- 1st provide a univaluated signal, chose one function from the provide_signals script
    1 chose a tag in the dictionary name_source
    2 if no tag exists, create it
    3 make sure the amplitude of the signal is renormalized to values vetween -2 to 2

- 2nd settle in the units of the signal by choosing a value for sampling_dt
    1 edit the 'unit' variable with a string bearing the physical unit
    2 make sure the sampling_dt corresponds to the unit selected!!!!
        correct the code to include if statements to make sure, for each dataset!

- 3rd do the fast fourier transfrom of the signal using the FFT(t, sig) function

- 4th select a scaling or frequency regime to compute the wavelets in these frequencies, these can be made
    1 by hand creating an arrray
    or 2 by selectring a range of frequencies
    or 3 by making it automatic using the last mode of the fourier transfer wc.scales_fourier(...)

- 4th compute the wavelets using one of the two possible methods by selecting the tag in wav_method, or both
    wav_method = 'pywt' -> wc.pywt_compute_wavelets(...)
    wav_method = 'niko' -> wc.niko_compute_wavelets(...)

- 5th select a name for the amplitude and phase of the wavelets to be stored or read
    1. select the name of the files
    2. make sure that the name corresponds to the dataset and analysis that you are about to use!!!

- 6th plot or compare the wavelets


IMPORTANT the most significant things in this script are 
    1 select the 'frequencies'
    2 make sure the 'units' and signal are 
'''



'''folder for data output'''
data_output_dir = './data/output/'  # '../data/' in case you are in the scripts folder

name_source = {
    'rain_india_manuel': './data/imput/s_allindiarain18712016.csv',
    'rossler_phase': './data/imput/rossler_phase_r12e27N03.dat',
    'ENSO_manuel': './data/imput/s_nino3418702020.csv',
    'ENSO_online': 'http://paos.colorado.edu/research/wavelets/wave_idl/sst_nino3.dat',
    'sin_signal': 'sin_gauss_'
}

'''tag for the frequencies bands and signal'''
#sig_tag = 'rossler_phase'
#sig_tag = 'rain_india_manuel'
sig_tag = 'ENSO_manuel'
#sig_tag = 'sin_signal'

freq_tag = 'lin'
#freq_tag = 'seed'
#freq_tag = 'log'


#unit = 'Hz'

def freq_generation(freq_tag = 'lin', step_period = 1 , min_period =1 , max_period =10, seeds =[1,10,100] ):
    '''frequency bands'''
    
    frequencies = []
    amplitudes = []
    if freq_tag == 'seed':

        seeds= [1/20., 1/100, 1/6, 1/3, 1/50,
                1/200]  # freq, they shall be below one
        #frequencies = 1/np.array([0.083, 0.1,0.5,0.9,1,2,4,5,6,7,8, 16, 32,64, 126])#
        amplitudes = [0.5, 1, 2, 1, 1, 2]
        frequencies = np.array(seeds)

    elif freq_tag == 'lin':

        step_period = 2 #months
        max_period = 85 #months
        min_period = 6 #months
        frequencies = 1./(np.arange(min_period, max_period, step_period))
        amplitudes = np.ones(len(frequencies))

    elif freq_tag == 'log':

        number_of_freq = 11
        max_period = 237
        min_period = 2
        period_ratio = np.log10(max_period)
        frequencies = 1/((min_period+np.arange(0, number_of_freq)**period_ratio))  
        amplitudes = np.ones(len(frequencies))

    if len(frequencies) < 1:
        raise Exception('you did not initialize well your tag CHECK! ')
    return frequencies, amplitudes


frequencies, amplitudes = freq_generation(freq_tag)
print ('ffff', frequencies[1:3], frequencies[-1], 'pppppp', 1/frequencies[1:3], 1/frequencies[-1])




'''depending on the the tag, load a signal and create time array'''
if sig_tag == 'rossler_phase':
    t, sig = ps.rossler_phase(name_source[sig_tag])
    unit = 'unit'

if sig_tag == 'ENSO_online':
    t, sig = ps.online_ENSO_34()
    unit = 'quar'

if 'ENSO_manuel' in sig_tag:
    t, sig = ps.read_ENSO_rain_manuel_files(name_source[sig_tag])
    sig = sig*20
    unit = 'mth'

if sig_tag == 'rain_india_manuel':
    t, sig = ps.read_ENSO_rain_manuel_files(name_source[sig_tag])
    unit = 'mth'

if sig_tag == 'sin_signal':
    '''compute signal'''
    gauss = False
    noise = False

    unit = 'Hz'

    sig_tag = name_source[sig_tag] + \
        str(gauss)[0] + 'noise_' + str(noise)[0]

    t = np.arange(600)
    t, sig = ps.create_signal(frequencies, amplitudes,
                              t, gauss=gauss, noise=noise)
    #sampling_dt = 1
    
    t, sig = ps.get_ave_values(t, sig, 3)

sampling_dt = t[1]-t[0]
#'''correct the signal time in spacific cases, depending on the units'''
#if unit == 'month' and 'manuel' in sig_tag:
#    t = np.arange(0, len(sig))

str_periods = freq_tag +'-p' + str(int(1/frequencies[0])) +\
     '-' + str(int(1/frequencies[-1])) + unit 

print('\n sampling period, last time and length of time', sampling_dt, t[-1], len(t), unit, '\n')


'''automatic frequecy scaling done by measuring the modes of the fourier transform'''
#par = wc.wavelets_scaling()
#bands_par = wc.wavelets_scaling(num_bands=6)#
#scales = wc.scales_fourier(wav.wavelet_kernel.fourier_factor(k0), par=par)

'''compute 1d fourier transformation'''
freq_spectrum, fft1d = wc.FFT(t, sig)#fft(sig)/len(t)
nyquist = int(len(fft1d))

'''compute wavelet decomposition for 2 different methods'''
wav_method = 'pywt'#'niko'#

if wav_method == 'pywt':
    kernel = 'cmor1.5-1.0'  # 'cmor'# #kind of wavelet kernel'gaussian'#
    waves, freq_bands = wc.pywt_compute_wavelets(sig, frequencies,
        kernel_name=kernel)
    
if wav_method == 'niko':
    kernel = 'morlet'
    waves, periods, freq_bands, cois = wc.niko_compute_wavelets(sig, frequencies,
        sampling_dt, kernel_name= kernel)


'''satore/read the amplitude and phase of the waveleets in/from numpy files'''
output_dir = data_output_dir + sig_tag + '_' + 'nSka_'+str(len(frequencies))+ '/'
name_files = sig_tag + '_' + wav_method + '_' + kernel + '_' + str_periods

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

wc.write_amplitude_phase_scale_wav(waves, 1.0 / frequencies, output_dir+name_files)
#amplitude, phase = wc.read_wavelets(name_files_pywt)

print(' \n saving amplitude ', output_dir+name_files + '_amp.npy')
print(' saving phase in ', output_dir+name_files + '_pha.npy')
print(' saving phase in ', output_dir+name_files + '_ska.npy\n')


'''call to create surrogates'''
n_surrogates = 111

for w, f  in zip(waves, freq_bands):
    surr_name = 'surr_circ_' + name_files
    #ident = 'f' + '{:04d}'.format(int(10000*f))
    surr_ident = '_p' + '{:04d}'.format(int(1/f))
    su.many_surrogates( output_dir + surr_name + surr_ident,  w,
                       min_shift=30, n_surrogates=n_surrogates)
    
print('\n surr stored withthis path name', output_dir + surr_name + surr_ident, '\n')


'''reconstruct the signal form the wavelets'''
#rec_signal_niko = wc.wav_reconstructed_signal(sig, waves, no_amp=False, individual_amp=True)


'''plot signals and wavelets'''
wp.plot_signal_phase_fft(t, sig, unit, freq_spectrum, frequencies, fft1d)
wp.plot_scalogram_fft_signal_together(
    t, sig, freq_bands,  waves, unit, waveletname=kernel)


#wp.plot_scalogram_fft_signal_together(
#    t, sig, freq_bands_pywt, waves_pywt, unit, waveletname=kernel_pywl)
#wp.plot_scalogram_fft_signal_together(
#    t, sig,freq_bands_pywt, waves_pywt, unit, waveletname=kernel_pywl)

#ax = wp.provide_axis()
#wp.plot_wavelet_scalogram(ax, t, freq_bands_niko, waves_niko, unit, waveletname = kernel_niko  )
#wp.plot_waves_amplitude_phase_WL('python wavelet', sig, rec_signal_pywt, waves_pywt, freq_bands_pywt )
#wp.plot_waves_amplitude_phase_WL('niko wavelet', sig, rec_signal_niko, waves_niko, freq_bands_niko)
#wp.plot_comparison_methods(waves_pywt[0], waves_niko[0], sig, rec_signal_pywt, rec_signal_niko)
wp.plot_all()

