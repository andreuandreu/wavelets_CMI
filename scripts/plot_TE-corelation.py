
from __future__ import division
from scipy.io import FortranFile
#from scipy.io import f90nml
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
import pickle
import matplotlib
import os
import matplotlib.ticker as ticker
 

#131072 samples of the first components of the coupled Roessler systems, for 100 different epsilons from 0 to almost 0.25.
#from zero eps to its maximum
#and compute CMI in both directions, using CMi as described in the paper, with 3-dim conditions
#obtained as x(t),x(t-5),x(t-10), where 5 and 10 are lags in number of samples
#python plot_TE-corelation.py ./correlations_TE_files.txt ./correlations_TE_files2.txt

def read_data(name):
    usecols = [0,1]
    print (name)
    #name = './correlations_TE/couplings_TE_VFb10-base2_1-5-10_122.txt'
    #x, y = np.loadtxt(name,  comments='#', delimiter=' ',
    #                  usecols=usecols, unpack=True)  # ,
    scales, TE = np.genfromtxt(name, delimiter=' ',
                               dtype="f8,f8", unpack=True, usecols=[0, 1])
    return scales, TE


def read_texts(fname ):

    f = open(fname,'r')

    names = []
    for line in f.readlines():
        names.append(line.replace('\n',''))

    f.close()

    return names


#ls = ['-', '-.','--',':', '-']



def plot_subfigures(ax, names):

    ls = ['-', ':','-', ':', '-', ':', '-', ':', '-', ':', '-', ':', '-', ':' ]
    ax.axvline(x=0.12)
    #ax.set_xlabel("Coupling strength $\epsilon$")
    #ax.set_ylabel("$I(a_0, b_{\u03C4}| b_0, b_5, b_{10} )$")
    
    for i, n in enumerate(names):
        x, y = read_data(n)
        #lb = n[-11:-4]
        ax.plot(x, y, label = n[-15:-4], ls = ls[i])
    ax.legend()
        #ax.legend(frameon=False, fontsize = 10 , loc = 'lower left')
   
    

fig = plt.figure('name')
ax = plt.subplot(1, 1, 1) 


names = read_texts(sys.argv[1]) 
ax = plt.subplot(1, 1, 1) 
plot_subfigures(ax, names)
    


#ax = fig.add_subplot(111)

plt.show()