
from __future__ import division
from cProfile import label
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

#python ./scripts/plot_TE-corelation.py 
def read_data(name):
    usecols = [0,1]

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

    name_files = sorted(fnmatch.filter(os.listdir(folder), root_name+'*') )

    size = len(name_files) 
    TEs_matrix = np.empty( shape=(size, size))
    for i, n in enumerate(name_files):
        #print('namTEma', i, n)
        TEs_matrix[i][:] = np.genfromtxt(folder+n, delimiter=' ',
                                   dtype=("f8"), unpack=True, usecols=[1])
        

    return TEs_matrix


def load_surrogateEntropies(folder, root_name):

    name_files = sorted( fnmatch.filter(os.listdir(folder), root_name+'*') )
    size = len(name_files)

    TEs_matrix = np.empty(shape=(size, 3, size))
    surr_TE_matrix = np.empty(shape=(size, size))
    for i, n in enumerate(name_files):
        #print(i, 'n', n)
        TEs_matrix[i][:] = np.genfromtxt(folder+"/"+n, delimiter=' ',
                                         dtype=("f8"), unpack=True, usecols=[2])
        surr_TE_matrix[i][:] = TEs_matrix[i][2][:]
    return TEs_matrix, surr_TE_matrix

def subtract_matrices(matA, matB):

    rA, cA = matA.shape
    subtraction = -matA + matB[0:rA, 0:cA]
    
    return subtraction

'''plot shit'''

def format_axes(ax, title, labels):

    ax.set_title(title, size=30)

    ax.tick_params(axis='both', which='major', labelsize=20)

    if 'ENSO' in  title or 'rain' in title:
        ax.set_xlabel(labels[0]+'[yr]', size=23)
        ax.set_ylabel(labels[1]+'[yr]', size=23)
        ax.xaxis.set_major_locator(MultipleLocator(12))
        ax.xaxis.set_major_formatter(
            FuncFormatter(lambda x, pos: int(x)/12))
        ax.xaxis.set_minor_locator(MultipleLocator(6))
        ax.yaxis.set_major_locator(MultipleLocator(12))
        ax.yaxis.set_major_formatter(
            FuncFormatter(lambda x, pos: int(x)/12))
        ax.yaxis.set_minor_locator(MultipleLocator(6))

    return ax

def plot_comoludogram(scales, phasPhas_TE, phasAmp_matrix):
    
    fig = plt.figure(figsize=(6, 6))

    xaxe, yaxe = np.meshgrid(scales, scales)

    gs = gridspec.GridSpec(1, 1)
    gs.update(left=0.05, right=0.95, hspace=0.3,
                top=0.95, bottom=0.05, wspace=0.15)
    axs = [gs[0, 0]]#, gs[0, 1]]
    toplot = [phasPhas_TE.T]#, phasAmp_matrix.T]
    #toplot = [phase_TE, phase_amp_TE]
    tits = ['PHASE-PHASE CAUSALITY']#, 'PHASE-AMP CAUSALITY']
    labs = ['PHASE', 'AMP']
    max_mat = max(map(max, phasAmp_matrix))
    for ax, cont, tit, lab in zip(axs, toplot, tits, labs):
        print('ttttt', type(ax))
        ax = plt.subplot(ax)
        cs = ax.contourf(xaxe, yaxe, cont, levels=np.arange(
            0.0, max_mat, 0.125), cmap=plt.cm.get_cmap("jet"), extend='max')
        # cs = ax.contourf(x, y, cont, levels = np.arange(4, 20, 0.125), cmap = plt.cm.get_cmap("jet"), extend = 'max')
        ax.tick_params(axis='both', which='major', labelsize=20)
        ax.set_title(tit, size=30)

        #ax = format_axes(ax, title, labs)
        
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
            plt.colorbar(cs)
            ax.grid()
            ax.set_ylabel("PERIOD %s [years]" % lab, size=23)
        
def plot_comoludogram_pixels(scales, matrix_TE, title, labels = ['pha', 'amp']):
    
    #fig, axs = plt.subplots(111)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    im = ax.imshow(matrix_TE, cmap="jet")
    plt.colorbar(im)
    ax = format_axes(ax, title, labels)



def plot_comoludogram_simple(scales, matrix_TE, title, labs = ['pha', 'amp']):

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    xaxe, yaxe = np.meshgrid(scales, scales)
    max_mat = max(map(max, matrix_TE))

    cs = ax.contourf(xaxe, yaxe, matrix_TE, levels=np.arange(
        0.0, max_mat, max_mat/10.), cmap=plt.cm.get_cmap("jet"), extend='max')
   
    ax = format_axes(ax, title, labs)
        
    plt.colorbar(cs)
    ax.grid()



    
#fig = plt.figure('name')
#ax = plt.subplot(1, 1, 1) 
#names = read_texts(sys.argv[1]) 
#ax = plt.subplot(1, 1, 1) 
#plot_subfigures(ax, names)
# sys.argv[1]


#folder = './data/output/corr_TE_ENSO_manuel_month_niko-cmor1_amp-amp/'
folder = './data/output/TE_ENSO_manuel_month_niko-p6-84mth/'
#root_name = 'Prob-est_VisFreq_b150'  # sys.argv[2]
root_name = 'Prob-est_VisFreq_b_bin-150_pha_pha_'#row-p6
#root_name = 'Prob-est_VisFreq_b_bin-150_pha_amp_'
surr_root = 'surr_circ_Prob-est_VisFreq_b_bin-150_SePha_Su-Pha_eDim-21111_'
#surr_root = 'surr_circ_Prob-est_VisFreq_b_bin-150_SeAmp_Su-Amp_eDim-21111_'
#surr_root = 'surr_circ_Prob-est_VisFreq_b_bin-150_amp_eDim-21111'
#surr_root = 'surr_circ_Prob-est_VisFreq_b_bin-150_pha_eDim-21111'
pha_or_amp = '_amp'



TEs_matrix = load_TransferEntropies(folder, root_name)
Surr_data, surr_TE_matrix = load_surrogateEntropies(folder, surr_root)
#subMat = subtract_matrices(surr_TE_matrix, TEs_matrix)

subMat = TEs_matrix  - surr_TE_matrix
step_period = 2 #months
max_period = 85 #months
min_period = 6 #months
scales = (np.arange(min_period, max_period, step_period))

#plot_comoludogram(scales, TEs_matrix, TEs_matrix)
title = 'ENSO Pha-Pha CMI'
#plot_comoludogram(scales, subMat, surr_TE_matrix)
labs = ['pha', 'pha']
#labs = ['pha', 'amp']
plot_comoludogram_pixels(scales, subMat, title, labs)


plot_comoludogram_simple(scales, subMat, title, labs)
#ax = fig.add_subplot(111)

plt.show()