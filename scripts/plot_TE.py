
import pandas as pd
import matplotlib.pyplot as plt
import sys
import numpy as np

from scipy import signal



def read_TransferEntropy_julia(name):
    ''' read formated file'''
    scales, TE = np.genfromtxt(name ,delimiter=' ' ,
                               dtype="f8,f8", unpack=True, usecols=[0, 1])

    return scales, TE

def plot_TE_rows(name, scales, TransferEntropy):

    fig = plt.figure('name')


    ax = plt.subplot(1, 1, 1)
    ax.plot(scales, TransferEntropy, label=name)
    ax.set_xscale('log')
    ax.legend(frameon=False, fontsize = 10 , loc = 'lower left')

name = sys.argv[1]
scales, TransferEntropy = np.genfromtxt(name, delimiter=' ',
                                        dtype="f8,f8", unpack=True, usecols=[0, 1])


plot_TE_rows(name[-15:-4], scales, TransferEntropy)
plt.show()