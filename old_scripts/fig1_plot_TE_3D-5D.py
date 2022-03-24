
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
    usecols = [0,1, 2]
    print (name)
    #name = './correlations_TE/couplings_TE_VFb10-base2_1-5-10_122.txt'
    x, y, sd = np.loadtxt(name,  comments='#', usecols = usecols, unpack = True )##, 

    return x, y, sd


def read_texts(fname ):

    f = open(fname,'r')

    names = []
    for line in f.readlines():
        names.append(line.replace('\n',''))

    f.close()

    return names



#ls = ['-', '-.','--',':', '-']


def plot_subfigures(i):

    ls = ['-', ':','-', ':', '-', ':', '-', ':', '-', ':', '-', ':', '-', ':' ]
    labels =  ['Kraskov x->y', 'Kraskov y->x'] #['VisitationFrequency x->y', 'VisitationFrequency y->x']#,KozachenkoLeonenko 
    ylab = ["$I(a_0, b_{\u03C4}| b_0 )$", "$I(a_0, b_{\u03C4}| b_0, b_5, b_{10} )$"]
    cl = ['r', 'b']
    ax = plt.subplot(2, 1, 1+i) 
    ax.axvline(x=0.12)
    ax.set_xlabel("Coupling strength $\epsilon$")
    ax.set_ylabel(ylab[i])
    
    for i, n in enumerate(names):
        x, y, sd = read_data(n)
        #lb = n[-11:-4]
        print(i, n)
        lb = labels[i]
        ax.plot(x, y, label = lb, ls = ls[i], color = cl[i])
        ax.fill_between(x,   y-sd, y+sd, color =  cl[i], alpha = 0.1)

    ax.legend()
    #ax.legend(frameon=False, fontsize = 10 , loc = 'lower left')
   
    

fig = plt.figure('name')

for (i,n) in enumerate(sys.argv[1:]):
    doc_files_name = n
    names = read_texts(doc_files_name) 
    plot_subfigures(i)
    


#ax = fig.add_subplot(111)

plt.show()