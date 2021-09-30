from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
import matplotlib
import os
import matplotlib.ticker as ticker
import math
from scipy.fftpack import *
from pactools import Comodulogram, REFERENCES
from pactools import simulate_pac


root = "/Users/andreu/Desktop/Dropbox/transfer_inormation_prague/"
dir = "TE_electrodes"
name = "/rossler_wavelets-r12e27N03.dat"
rosslerData = np.loadtxt(root+dir+name, dtype= float)

def plot_basic():

    plt.plot(rosslerData[0:100])
    plt.plot(ht[0:100])


if __name__ == '__main__':
    N = len(rosslerData)
    f = 1
    dt = 1.0 / N
    
    ht = hilbert(rosslerData)

