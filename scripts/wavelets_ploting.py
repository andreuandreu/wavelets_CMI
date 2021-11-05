from numpy.core.numeric import base_repr
from numpy.testing._private.utils import nulp_diff
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


import numpy as np
from numpy.fft import fft, ifft

import provide_signals as css

'''
 plots
        plot_signal_phase_fft(time, signal): signal + phase transform + FFt in frequency and 
        wavenumber domain
        plot_waves_amplitude_phase_WL(): reconstruction of the signal and each wavelet 
        (phase and amplitude) corresponding to the provided frequency
        plot_comparison_methods(): to compare the wavelets and reconstructions of two 
        different wavelet methods ->in development<-   
        plot_scalogram_fft_signal_together():
        plot_scalogram():

'''


def plot_signal_plus_average(ax, time, signal, unit, smooth_len= 5):

    # plotting the signal
    ax.set_xlabel('Time [' + unit + ']')
    ax.set_ylabel('Amplitude')
    ax.set_title("Signal + Time Average")
    ax.set_xlim([time[0], time[-1]])
    ax.plot(time, signal, color='green')

    #plotting the average of the signal
    time_ave, signal_ave = css.get_ave_values(time, signal, smooth_len)
    ax.plot(time_ave, signal_ave, label='time average (n={})'.format(smooth_len))
    ax.legend(loc='upper right', fontsize='small', frameon=False)


def plot_signal_phase_fft(time, signal, unit, freq_spectrum, frequencies, fft1d):

    plt.rcParams['figure.figsize'] = [12, 7]
    fig, ax = plt.subplots(3, 1)

    # plotting the signal
    plot_signal_plus_average(ax[0], time, signal, unit)
  
    # plotting the phase spectrum of the signal
    ax[1].set_title("Phase Spectrum of the Signal")
    ax[1].phase_spectrum(signal, color='green')

    # plotting the fft of the signal in frequency domain
    ax[2].set_title("fft Spectrum of the Signal")
    #ax[2].loglog(fft1d[0:nyquist], 'r')
    ax[2].plot(freq_spectrum, 2*fft1d/len(time))
    ax[2].scatter(frequencies, np.ones(len(frequencies)))
    ax[2].set_yscale('log')
    ax[2].set_xscale('log')

    # plotting the fft of the signal in scale/length and frequency domain and log scale
    ax2 = ax[2].twiny()
    ax2.plot(1./freq_spectrum, fft1d/len(time), c='r')
    #ax2.scatter(1./np.array(freq_bands_niko), np.ones( len(freq_bands_niko) ), c='r')
    #ax2.semilogx(len(t)/freq_spectrum**2, fft1d/len(t), c='r')
    ax2.set_xlabel('wavenum', color='r')
    #ax[2].loglog(frequencies, 'b')
    plt.tight_layout()


def plot_fft(ax, freq_spectrum, fft1d, plot_direction='vertical', yticks=None, ylim=None):

    if plot_direction == 'vertical':

        #fft with scales (periods/dimension) in log scale
        scales = 1./freq_spectrum
        scales_log = np.log2(scales)
        ax.scatter([np.mean(fft1d)-0.5] * len(freq_spectrum), np.log2(1/freq_spectrum))
        ax.plot(fft1d/np.mean(fft1d), scales_log,
                'r-', label='Fourier Transform')
        ax.legend(loc='upper right', fontsize='small', frameon=False)
        
        #set log scale
        ax.set_yticks(np.log2(yticks))
        ax.set_yticklabels(yticks)
        ax.invert_yaxis()
        ax.set_ylim(ylim[0], ylim[-1])

        #fft with frequencies
        ax2 = ax.twinx()
        ax2.plot(fft1d/np.mean(fft1d), freq_spectrum, 'g--')
        #ax2.set_yscale('linear')
        ax2.set_ylabel('Frequency', color='g')


def plot_waves_amplitude_phase_WL(title, sig, rec_sig, waves, frequencies):

    fig, axs = plt.subplots(nrows=3, ncols=len(
        waves),  sharey="row", figsize=(15, 10))  # sharex=True,
    fig.suptitle(title, fontsize=20)

    axs[0, 0].set_ylabel("signal")
    axs[1, 0].set_ylabel("amplitude")
    axs[2, 0].set_ylabel("phase")

    for i in range(len(waves)):
        axs[0, i].set_title(' F  {0:.2f}  P {1:.2f}'.format(
            frequencies[i], 1./frequencies[i]))
        axs[0, i].plot(sig[0: int(len(waves[i])/(i+1))], label='original')
        axs[0, i].plot(rec_sig[0: int(len(waves[i])/(i+1))],
                       label='reconstructed')
        axs[1, i].plot(np.real(waves[i, :]))
        axs[1, i].plot(np.abs(waves[i, :]))
        axs[2, i].plot(np.angle(waves[i, :]))
    axs[0, 0].legend(loc='upper right', fontsize='small', frameon=False)


def plot_comparison_methods(wav1, wav2, sig, rec_signal_1, rec_signal_2):

    fig2, ax = plt.subplots(4, 1)
    # plotting the amplitude
    ax[0].set_xlim(20, len(wav1)-20)
    ymax = max(wav1[0][20:-20]) + 0.3*max(wav2[0][20:-20])
    #ax[0].set_ylim( -ymax, ymax   )
    ax[0].set_title("real wavelets")
    ax[0].plot(wav1.real)
    #ax[0].plot(wav1[0].imag)
    #ax[0].plot(wav2[0].imag)
    ax[0].plot(wav2.real)

    # plotting the phase
    ax[1].set_xlim(20, len(wav1)-20)
    ax[1].set_title(" phase whavelets")

    ax[1].plot(np.angle(wav1))
    ax[1].plot(np.angle(wav2))

    # plotting the amp
    ax[2].set_xlim(20, len(wav1)-20)
    ax[2].set_title(" amp whavelets")
    ax[2].plot(np.abs(wav1))
    ax[2].plot(np.abs(wav2))

    # plotting the reconstruction
    ax[3].set_xlim(20, 180)
    ax[3].plot(sig, color='green')
    ax[3].plot(rec_signal_1, color='blue')
    ax[3].plot(rec_signal_2, color='red')


def plot_scalogram_fft_signal_together(time, signal, freq, wavelets, unit, waveletname):

    #prepare plot
    fig = plt.figure(figsize=(12, 12))
    spec = gridspec.GridSpec(ncols=6, nrows=6)
    top_ax = fig.add_subplot(spec[0, 0:5])
    bottom_left_ax = fig.add_subplot(spec[1:, 0:5])
    bottom_right_ax = fig.add_subplot(spec[1:, 5])

    #plot signal
    plot_signal_plus_average(top_ax, time, signal, unit)

    #plot scalogram
    yticks, ylim, im = plot_wavelet_scalogram(
        time, freq, wavelets, unit, waveletname, bottom_left_ax)
    #position colorbar
    cbar_ax = fig.add_axes([0.9, 0.85, 0.01, 0.1])
    fig.colorbar(im, cax=cbar_ax, orientation="vertical")

    #plot fft
    plot_fft(bottom_right_ax, plot_direction='vertical',
             yticks=yticks, ylim=ylim)
    bottom_right_ax.set_ylabel('Period [' + unit + ']',  fontsize=14)
    plt.tight_layout()


def plot_wavelet_scalogram(time, freq, wavelets, unit, waveletname='cmor', ax=plt.subplots(figsize=(15, 10))):

    #prepare data
    power = (abs(wavelets)) ** 2
    period = 1.0 / (freq)  # 1/freq
    #prepare contours levels
    # [0.0625, 0.25,  1, 8]##[8, 16000, 32000, 128000]#
    levels = [0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8]
    contourlevels = np.log2(levels)
    cmap = plt.cm.seismic
    im = ax.contourf(time, np.log2(period), np.log2(power),
                     contourlevels, extend='both', cmap=cmap)

    #prepare axis tags
    ax.set_title('Wavelet '+waveletname +
                 ' Transform/time Power Spectrum)', fontsize=20)
    ax.set_ylabel('Scale [' + unit + ']', fontsize=18)
    ax.set_xlabel('Time [' + unit + ']', fontsize=18)

    #arrange y axis scales
    yticks = 2**np.arange(np.ceil(np.log2(period.min())),
                          np.ceil(np.log2(period.max())))
    ax.set_yticks(np.log2(yticks))
    ax.set_yticklabels(yticks)
    ax.invert_yaxis()
    ylim = ax.get_ylim()
    ax.set_ylim(ylim[0], ylim[-1])

    return yticks, ylim, im


fig, ax = plt.subplots(1)
ax.plot(freq_spectrum)
ax.plot(np.arange(0, len(periods_niko))*len(freq_spectrum) /
        len(periods_niko), 1/periods_niko[::-1])


'''plot signals and wavelets'''
plot_signal_phase_fft(t, sig, freq_spectrum, frequencies, fft1d)
plot_scalogram_fft_signal_together(
    t, sig,  1/periods_niko, waves_niko, unit, waveletname=kernel_niko)
#plot_scalogram_fft_signal_together(
#    t, sig, freq_bands_pywt, waves_pywt, unit, waveletname=kernel_pywl)


plot_signal_phase_fft(t, sig, freq_spectrum, frequencies, fft1d)
plot_scalogram_fft_signal_together(t, sig, 1/periods_niko, waves_niko, unit, waveletname = kernel_pywl)

#plot_scalogram_fft_signal_together(
#    t, sig,freq_bands_pywt, waves_pywt, unit, waveletname=kernel_pywl)

#plot_wavelet_scalogram(t, freq_bands_niko, waves_niko, waveletname = kernel_niko  )
#plot_wavelet_scalogram(t, freq_bands_pywt, waves_pywt, waveletname = kernel_pywl  )
#plot_waves_amplitude_phase_WL('python wavelet', sig, rec_signal_pywt, waves_pywt, freq_bands_pywt )
#plot_waves_amplitude_phase_WL('niko wavelet', sig, rec_signal_niko, waves_niko, freq_bands_niko)
#plot_comparison_methods()
plt.show()