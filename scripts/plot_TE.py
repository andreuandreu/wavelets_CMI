
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

def plot_TE_rows(name,  TransferEntropy):

    fig = plt.figure('name')


    ax = plt.subplot(1, 1, 1)
    t = (TransferEntropy -np.mean(TransferEntropy) ) / np.std(TransferEntropy)
    ax.plot(t, label=name)
    #ax.set_xscale('log')
    ax.legend(frameon=False, fontsize = 10 , loc = 'lower left')

names = sys.argv[:]

for n in names:

    #scales, TransferEntropy = np.genfromtxt(n, delimiter=' ',
    #                                    dtype="f8,f8", unpack=True, usecols=[0, 1])
    TransferEntropy = np.genfromtxt(n, delimiter=' ',
                                            dtype="f8", unpack=True, usecols=[0])

    plot_TE_rows(n[-15:-4], TransferEntropy)
plt.show()