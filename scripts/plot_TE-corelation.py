
from __future__ import division
from cProfile import label
from logging import exception
from turtle import color
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
import create_config as cc

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

def load_TransferEntropies(folder, data_root):

    

    name_files = sorted(fnmatch.filter(os.listdir(folder), data_root+'*') )
    

    print(name_files)
    size = len(name_files) 
    TEs_matrix = np.empty( shape=(size, size))
    for i, n in enumerate(name_files):
        print('namTEma', i, n)
        TEs_matrix[i][:] = np.genfromtxt(folder+n, delimiter=' ',
                                   dtype=("f8"), unpack=True, usecols=[1])
        

    return TEs_matrix


    


def load_surrogateEntropies(folder, data_root):

    name_files = sorted( fnmatch.filter(os.listdir(folder), data_root+'*') )

    size = len(name_files)

    if size < 1:
        print('nothing here ', folder+data_root)
        print('ERRROR, bad path to names')
        raise SystemExit
        

    TEs_matrix = np.empty(shape=(size, 2, size))
    surr_TE_matrix = np.empty(shape=(size, size))
    surr_sdTE_matrix = np.empty(shape=(size, size))
    all_surr = np.array([])
    for i, n in enumerate(name_files):
        print(i, 'n', n)
        TEs_matrix[i][:][:] = np.genfromtxt(folder+"/"+n, delimiter=' ',
                                         dtype=("f8", "f8"), unpack=True, usecols=[2, 3])
        surr_TE_matrix[i][:] = TEs_matrix[i][0][:]
        surr_sdTE_matrix[i][:] = TEs_matrix[i][1][:]

        np.append( all_surr, np.load(folder + "set-" + n[:-4] + '.npy' ))
    
    return surr_TE_matrix, surr_sdTE_matrix

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


def plot_projected_TE_vs_surr(scales, data_arr, data_var, surr_arr, surr_var ):

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    ax.plot(scales, data_arr, color = 'r' )
    ax.errorbar(scales, data_arr, yerr = data_var, color = 'r')

    ax.plot(scales, surr_arr, color = '0.8' )
    ax.errorbar(scales, surr_arr, yerr = surr_var, color = '0.8')



def plot_projected_diff(scales, arr, var ):

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    
    ax.plot(scales, arr, color = 'r', alpha = 0.5 )
    ax.errorbar(scales, arr, yerr = var, color = 'r', alpha = 0.5)




def plot_data_sd(scales, dat_arrays, colors, labels ):

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    
    ax.axhline(y=0, color='black', linestyle='-')
    ax.set_ylabel('TE-TEsurr')
    ax.set_xlabel('period')
    shift = [0, 0.4]

    for e, c, l, s in zip(dat_arrays, colors, labels, shift):
        ax.plot(scales+s, e[0], alpha = 0.5 , color = c, label = l)
        ax.errorbar(scales+s, e[0], yerr =e[1], alpha = 0.5, color = c )
    ax.legend(frameon = False)



def projection(matrix):

    values_array = np.empty(shape=(matrix.shape[0]))
    std_array = np.empty(shape=(matrix.shape[0]))
    for i, c in enumerate(matrix):
         values_array[i] = np.mean(c)
         std_array[i] = np.std(c)
    
    return values_array, std_array



def matrixflip( m, d = 'h'):
    myl = np.array(m)
    if d=='v': 
        return np.flip(myl, axis=0)
    elif d=='h':
        return np.flip(myl, axis=1)



tag = 'sin'
#to_plot = ['pha', 'amp']
to_plot = ['pha', 'pha']
#wavelet = 'niko'

#name_config_file = './confs/config_embeding_char_ENSO_pywt.ini'
name_config_file = './confs/config_embeding_char.ini'
conf = cc.load_config(name_config_file )
folder = conf['folders']['data_folder']+ conf['folders']['export_folder'] 
root_name =  conf['prob_est']['name_tag']+"_bin-" +  conf['emb_par']['bins']


if to_plot[0] == 'pha' and to_plot[1] == 'pha':
    '''pha-pha '''
    data_root = root_name + '_pha_pha_row-'
    surr_root = root_name + '_SePha_Su-Pha_eDim-21111_'
    labels = ['pha', 'pha']

if to_plot[0] == 'pha' and to_plot[1] == 'amp':
    '''pha-amp '''
    data_root = root_name + '_pha_amp_row-'
    surr_root = root_name + '_SePha_Su-Amp_eDim-21111_'
    labels = ['pha', 'amp']

labels = to_plot
title = 'CMI ' + conf['names']['in_data_tag'][0:4] + ' ' + labels[0] + '-' + labels[1]

folder_TEdata = folder + "Dat_files/"
TEs_matrix = load_TransferEntropies(folder_TEdata, data_root)

folder_surrogates = folder + "Su-"+conf['surrogates']['surr_kind'] + "_files/"
surr_TE_matrix, surr_sdTE_matrix = load_surrogateEntropies(folder_surrogates, surr_root)

subMat = TEs_matrix  - surr_TE_matrix

data_arr, data_var = projection(TEs_matrix)
surr_arr, surr_var = projection(surr_TE_matrix)

dif_arr, dif_var = projection(subMat.T)



max_period = 85 #months
min_period = 6 #months

step_period = int( (max_period - min_period) / np.shape(subMat)[0] + 0.5 )
scales = (np.arange(min_period, max_period, step_period))

plot_comoludogram_pixels(scales, matrixflip(subMat, 'v'), title, labels)

plot_comoludogram_simple(scales, subMat, title, labels)

#plot_projected_TE_vs_surr(scales, data_arr, data_var, surr_arr, surr_var )

data_arrays = [projection(subMat), projection(subMat.T)]

colors = ['r', 'b']
labels = ['over x', 'over y']

#plot_projected_diff(scales, dif_arr, dif_var )
plot_data_sd(scales, data_arrays, colors, labels )


plt.show()