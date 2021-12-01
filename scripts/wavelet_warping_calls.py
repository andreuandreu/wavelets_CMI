import numpy as np
import provide_signals as ps
import wavelets_computation as wc
import wavelets_ploting as wp
import configparser

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
data = './data/output/'  # '../data/' in case you are in the scripts folder



gauss = False
noise = False
name_source = {
    'rain_india_manuel': './data/imput/s_allindiarain18712016.csv',
    'rossler_phase': './data/imput/rossler_phase_r12e27N03.dat',
    'ENSO_manuel': './data/imput/s_nino3418702020.csv',
    'ENSO_online': 'http://paos.colorado.edu/research/wavelets/wave_idl/sst_nino3.dat',
    'synthetic': 'sin_gauss_' + str(gauss)[0] + 'noise_' + str(noise)[0]
}

'''compute signal'''
seed_freq = [1/20., 1/100, 1/6, 1/3, 1/50, 1/200]  # freq, they shall be below one
amplitudes = [0.5, 1, 2, 1, 1, 2]
#sampling_dt = 1
t = np.arange(600)

t, sig = ps.create_signal(seed_freq, amplitudes, t, gauss=gauss, noise=noise)
#sig_tag = 'synthetic'
sig_tag = 'rossler_phase'
t, sig = ps.rossler_phase(name_source['rossler_phase'])

'''call for the time and signal'''
#sig_tag = 'rain_india_manuel'
#sig_tag = 'ENSO_manuel'
#sig_tag = 'sin_signal'
#t, sig = ps.read_ENSO_rain_manuel_files(name_source[sig_tag])
#t, sig = ps.online_ENSO_34()
if sig_tag == 'ENSO_manuel':
    sig = sig*20
#t, sig = ps.get_ave_values(t, ig, 3)



'''characteristics of the signal and the processing'''
#unit = 'month'
unit = 'Hz'

'''correct the signal time in spacific cases, depending on the units'''
if unit == 'month' and 'manuel' in sig_tag:
    t = np.arange(0, len(sig))

sampling_dt = t[1]-t[0]
print('\n sampling period, last time and length of time', sampling_dt, t[-1], len(t), '\n')
#frequencies = 1/np.array([0.083, 0.1,0.5,0.9,1,2,4,5,6,7,8, 16, 32,64, 126])#
if sig_tag == 'synthetic':
    frequencies = np.array(seed_freq)
else:
    frequencies = 1/((np.arange(1, 72)*sampling_dt))  #

'''automatic frequecy scaling done by measuring the modes of the fourier transform'''
#par = wc.wavelets_scaling()
#scales = wc.scales_fourier(wav.wavelet_kernel.fourier_factor(k0), par=par)

'''compute 1d fourier transformation'''
freq_spectrum, fft1d = wc.FFT(t, sig)#fft(sig)/len(t)
nyquist = int(len(fft1d))


'''compute wavelet decomposition for 2 different methods'''
kernel_pywl = 'cmor1.5-1.0'  # 'cmor'# #kind of wavelet kernel'gaussian'#
kernel_niko = 'morlet'
wav_method = 'niko'#pywt
bands_par = wc.wavelets_scaling(num_bands=6)

if wav_method == 'pywt':
    waves, freq_bands = wc.pywt_compute_wavelets(sig, frequencies,
        kernel_name=kernel_pywl)
if wav_method == 'niko':
    waves, periods, freq_bands, cois = wc.niko_compute_wavelets( sig, frequencies,
        sampling_dt, kernel_name =  kernel_niko)


'''satore/read the amplitude and phase of the waveleets in/from numpy files'''
name_files = data + sig_tag + '_' + 'Nska_'+str(len(frequencies))+ unit + '_' + wav_method + '_' + kernel_pywl
#wc.write_amplitude_phase_wav(waves_pywt, name_files_pywt)
wc.write_amplitude_phase_scale_wav(waves, 1.0 / frequencies, name_files)
#amplitude, phase = wc.read_wavelets(name_files_pywt)
print("\n npy files stored in ", name_files, '\n')

'''reconstruct the signal form the wavelets'''
rec_signal_niko = wc.wav_reconstructed_signal(sig, waves, no_amp=False, individual_amp=True)

'''plot signals and wavelets'''
wp.plot_signal_phase_fft(t, sig, unit, freq_spectrum, frequencies, fft1d)
wp.plot_scalogram_fft_signal_together(
    t, sig,  1/periods,  waves, unit, waveletname=kernel_niko)
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
