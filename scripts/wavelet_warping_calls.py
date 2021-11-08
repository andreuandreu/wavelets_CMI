
import wavelet_analysis as wa
import numpy as np
from numpy.fft import fft, ifft
import pywt
from scipy import signal
import provide_signals as ps
import wavelets_computation as wc


'''

- 1st provide a univaluated signal, chose one function from the provide_signals script

- 2nd settle in the units of the signal by chosing a value for sampling_dt

- 3rd do the fast fourier transfrom of the signal using the FFT(t, sig) function

'''



'''folder for data output'''
data = './data/output/'  # '../data/' in case you are in the scripts folder

'''compute signal
seed_freq = [1/20., 1/100, 1/6]  # freq, they shall be below one
amplitudes = [0.5, 1, 2]
#sampling_dt = 1
t = np.arange(600)

t, sig = ps.create_signal(seed_freq, amplitudes, t, gauss=gauss, noise=noise)
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
t, sig = ps.read_ENSO_rain_manuel_files(name_source[sig_tag])
#t, sig = ps.online_ENSO_34()

if sig_tag == 'ENSO_manuel':
    sig = sig*20

time = np.arange(0, len(sig))
#t, sig = ps.get_ave_values(t, sig, 3)

'''characteristics of the signal and the processing'''
unit = 'month'
sampling_dt = 1# t[1]-t[0]
print('\n sampling period, last time and length of time', sampling_dt, time[-1], len(time), '\n')
#frequencies = 1/np.array([0.083, 0.1,0.5,0.9,1,2,4,5,6,7,8, 16, 32,64, 126])#
frequencies = 1/((np.arange(1, 256)*sampling_dt))  #

#par = wc.wavelets_scaling()
#scales = wc.scales_fourier(wav.wavelet_kernel.fourier_factor(k0), par=par)
'''compute 1d fourier transformation'''
freq_spectrum, fft1d = wc.FFT(time, sig)#fft(sig)/len(t)
nyquist = int(len(fft1d))


'''compute wavelet decomposition for 2 different methods'''
kernel_pywl = 'cmor1.5-1.0'  # 'cmor'# #kind of wavelet kernel'gaussian'#
kernel_niko = 'morlet'
bands_par = wc.wavelets_scaling(num_bands=6)
waves_pywt, freq_bands_pywt = wc.pywt_compute_wavelets(sig, frequencies,
    kernel_name=kernel_pywl)
waves_niko, periods_niko, freq_bands_niko, cois = wc.niko_compute_wavelets( sig, frequencies,
    kernel_name =  kernel_niko)


'''satore the amplitude and phase of the waveleets in numpy files'''
name_files_pywt = data + sig_tag +'_'+ unit +'_wavelet_vecors_pywt_'+ kernel_pywl 
wc.write_amplitude_phase_wav(waves_pywt, name_files_pywt)


'''reconstruct the signal form the wavelets'''
rec_signal_pywt = wc.wav_reconstructed_signal(sig, waves_pywt, no_amp=False, individual_amp=True)
rec_signal_niko = wc.wav_reconstructed_signal(sig, waves_niko, no_amp=False, individual_amp=True)