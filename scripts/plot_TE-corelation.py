
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
from matplotlib.ticker import MultipleLocator, FuncFormatter
import matplotlib.gridspec as gridspec
import fnmatch

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

def load_TransferEntropies(folder, root_name):

    name_files = fnmatch.filter(os.listdir(folder), root_name+'*')
    

    size = len(name_files) -1

    TEs_matrix = np.empty( shape=(size, size))
    for i, n in enumerate(name_files[:-1]):

        TEs_matrix[i][:] = np.genfromtxt(folder+n, delimiter=' ',
                                   dtype="f8", unpack=True, usecols=[0])
    return TEs_matrix



def plot_comoludogram(scales, phase_TE, phase_amp_TE):
    
    fig = plt.figure(figsize=(15, 7.5))

    xaxe, yaxe = np.meshgrid(scales, scales)

    gs = gridspec.GridSpec(1, 2)
    gs.update(left=0.05, right=0.95, hspace=0.3,
                top=0.95, bottom=0.05, wspace=0.15)
    axs = [gs[0, 0], gs[0, 1]]
    toplot = [phase_TE.T, phase_amp_TE.T]
    #toplot = [phase_TE, phase_amp_TE]
    tits = ['PHASE-PHASE CAUSALITY', 'PHASE-AMP CAUSALITY']
    labs = ['PHASE', 'AMP']

    for ax, cont, tit, lab in zip(axs, toplot, tits, labs):
        ax = plt.subplot(ax)
        cs = ax.contourf(xaxe, yaxe, cont, levels=np.arange(
            0.99, 1, 0.00125), cmap=plt.cm.get_cmap("jet"), extend='max')
        # cs = ax.contourf(x, y, cont, levels = np.arange(4, 20, 0.125), cmap = plt.cm.get_cmap("jet"), extend = 'max')
        ax.tick_params(axis='both', which='major', labelsize=20)
        ax.set_title(tit, size=30)

        if 'ENSO' in folder or 'rain' in folder:
            ax.xaxis.set_major_locator(MultipleLocator(12))
            ax.xaxis.set_major_formatter(
                FuncFormatter(lambda x, pos: int(x)/12))
            ax.xaxis.set_minor_locator(MultipleLocator(6))
            ax.yaxis.set_major_locator(MultipleLocator(12))
            ax.yaxis.set_major_formatter(
                FuncFormatter(lambda x, pos: int(x)/12))
            ax.yaxis.set_minor_locator(MultipleLocator(6))
            ax.set_xlabel("PERIOD PHASE [years]", size=23)
            # plt.colorbar(cs)
            ax.grid()
            ax.set_ylabel("PERIOD %s [years]" % lab, size=23)
    #plt.savefig('plots/nino34-CESMhigh.eps', bbox_inches="tight")

    

#fig = plt.figure('name')
#ax = plt.subplot(1, 1, 1) 
#names = read_texts(sys.argv[1]) 
#ax = plt.subplot(1, 1, 1) 
#plot_subfigures(ax, names)

folder = sys.argv[1]
root_name = sys.argv[2]

TEs_matrix = load_TransferEntropies(folder, root_name)
scales = np.arange(len(TEs_matrix))
plot_comoludogram(scales, TEs_matrix, TEs_matrix)

#ax = fig.add_subplot(111)

plt.show()