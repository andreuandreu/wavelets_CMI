import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
from glob import glob
from matplotlib import cm
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.ndimage.filters import gaussian_filter
import scipy.ndimage
from scipy.interpolate import griddata
import math
import read_matlab_fig as rmf



def read_milan_xyy_4col(names):

    data = np.empty((4, 2, 50))

    for i, n in enumerate(names):
        d = np.genfromtxt(n, dtype=("f8", "f8" , "f8" , "f8"), unpack=True, usecols=[1, 2, 3, 4])
        data[2*i] = d[0:2]
        data[2*i+1] = d[2:4]
    
    return data

def read_milan_xyy_2col(name):

    data = np.genfromtxt(name, dtype=("f8", "f8" ), unpack=True, usecols=[1, 2])
    return data

def read_plot_milan_rat_lags(data_names, axes):

    xpos = [26,75,75,26]
    ymin = 0
    ymax = 1

    for key, n, xpo in zip( axes, data_names, xpos):
        print ('rat', n, key)
 
        labels = ['', '']
        data = read_milan_xyy_2col(n)
        
        if key == 'upper left':
            labels = [r'HF amp $\rightarrow$ LF pha', r'LF pha $\rightarrow$ HF amp']

        if key == 'upper right':
                labels = ['HF amp leads', 'LF pha leads']

        plot_milan_lags(axes[key], data, labels)

        if 'right' in key:
            axes[key].axvline( xpo, ymin, ymax)
            axes[key].set_yticks([])
        if 'upper' in key:
            axes[key].set_xticks([])

def read_plot_milan_ross_lags(data, axes):

    xpos = [ 13, 10,  10, 13]

    for key, d, pos in zip(axes, data, xpos):
        
        labels = ['', '']
        print ('ross',  key)
        
        axes[key].set_xlim(0, 33)

        if 'upper' in key:
            axes[key].set_xticks([])

        if key == 'upper left':
            labels = [r'y $\rightarrow$ x', r'x $\rightarrow$ y']

        if key == 'lower right':
            labels = ['y leads', 'x leads']

        if 'left'in key:
            plot_milan_lags(axes[key], d, labels)
            

        if 'right' in key:
            plot_milan_lags_ross(axes[key], d, pos, labels)
            axes[key].set_yticks([])

   
def filter_data( Z, cut = 2.0, filter = 3,  percentage = 0.66):

    #thresholdZ = max(Z) - (max(Z)-min(Z)) * percentage
    filtered_z = scipy.ndimage.zoom(Z, filter)#griddata(Z,  method='nearest')#gaussian_filter(Z, sigma=filter)#
    
    masked_array = np.ma.masked_where( (filtered_z > -cut) & (filtered_z < cut) , filtered_z)
    
    return masked_array           


def plot_milan_lags_ross(ax, data, xpos, labels, ymin = 0, ymax = 1):

    ax.set_xlim(0, 33)
    ax.axvline( xpos, ymin, ymax)
    ax.plot( (data[0]) / (np.max(data[1])-np.min(data[1]) ), color = 'g', label = labels[0])
    ax.plot( (data[1]) / (np.max(data[1])-np.min(data[1]) ), color = 'r', label = labels[1])

    ax.legend(frameon=False)

def plot_milan_lags(ax, data, labels):

    ax.plot( (data[0]) / (np.max(data[1])-np.min(data[1]) ), color = 'g',  label = labels[0])
    ax.plot( (data[1]) / (np.max(data[1])-np.min(data[1]) ), color = 'r',  label = labels[1] )

    ax.legend(frameon=False)
    

def plot_milan_xyz(ax, dataX, dataY, dataZ):

    lenTot = len(dataX)
    lenX = 71
    lenY = int(lenTot/lenX)
    zAux = dataZ.reshape(lenY, lenX)
    
    filter = 12
    filteredZ = filter_data(zAux, cut = cut, filter = filter)

    axeXpoints = []
    for i in range(lenY):
        axeXpoints.append(dataX[i*lenX])

    lin = np.arange(lenY)
    Xinterp = scipy.interpolate.interp1d(lin, axeXpoints)
    lin2 = np.linspace(0, lenY-1 , num=lenY*filter, endpoint=True ) 
    axeXpoints = Xinterp(lin2)

    lin = np.arange(lenX)
    Yinterp = scipy.interpolate.interp1d(lin, dataY[0:lenX])
    lin2 = np.linspace(0, lenX-1, num=lenX*filter, endpoint=True )
    axeYpoints = Yinterp(lin2)
    
    cmap = plt.get_cmap('jet')#''YlGn'
    cmap.set_bad(color='white')
    #im = ax.imshow(  np.flipud(zAux.T), aspect = 'auto', cmap = cmap  ) #axeXpoints, dataY[0:lenX],
    #interpolation='gaussian' # extent=(np.amin(dataX), np.amax(dataX), np.amin(dataY), np.amax(dataY)),
    im = ax.contourf( axeXpoints, axeYpoints ,  filteredZ.T,
                      extend='both', cmap=cmap)#contourlevels, 

    return(im)

def prepare_figue_lag(fig, axes):


    fig.text(0.25, 0.9, 'CMI', va='center', fontsize=13)
    fig.text(0.65, 0.9, 'MI', va='center', fontsize=13)

    fig.text(0.5, 0.04, 'Lag[su]', va='center', fontsize=11)
    fig.text(0.022, 0.5, 'I[nu]', va='center', rotation='vertical', fontsize=11)  

    #fig.text(0.5, 0.04, '8Hz, HF amp 80 Hz', va='center', fontsize=11)
    #fig.text(0.5, 0.04, '8Hz, HF amp 100 Hz', va='center', fontsize=11)
    #IM-IC LF 8Hz, HF 100 HZ
    #fig.text(0.02, 0.2, 'PP-IC', va='center', rotation='vertical', fontsize=13)
    #fig.text(0.02, 0.5, 'Im-IC', va='center', rotation='vertical', fontsize=13)


#x10 v1 s1 phAA 12k1LX3ct1 kL4 shi1 su30 v1, 2, 3 variables (Sch-IC, Im-IC, PP-IC) 
#x10 v1 s1 AAph 21k1LX3ct1 kL4 shi1 su30


#filename = "../../Desktop/Dropbox/transfer_inormation_prague/plots/fig1_lag/s10v2a8f50TPL3.xyy"
#data = read_milan_xyy_2col(filename)
##plot_milan(data)

filename = "../../Desktop/Dropbox/transfer_inormation_prague/plots/fig1_lag/r11e31phTPcd4mi2q4.xyy"
#data = read_milan_xyy(filename)
#plot_milan(data)



'''FIGURE Rat data lag'''
pth1 ="../../Desktop/Dropbox/transfer_inormation_prague/plots/fig1_lag/"
root = "s10*.xyy"
data_names = glob(pth1+root)
mossaic_keys = [['upper left', 'upper right'],
                ['lower left', 'lower right']]
fig, axes = plt.subplot_mosaic(mossaic_keys,  constrained_layout=True, 
                               gridspec_kw={'hspace': 0, 'wspace': 0})#figsize=(5.5, 3.5)
fig.subplots_adjust(hspace=0)
read_plot_milan_rat_lags(data_names, axes)
prepare_figue_lag(fig, axes)
plt.savefig(pth1+"fig1_rat_lags.svg")

'''FIGURE Ross data lag'''
root = "r11*.xyy"
data_names = glob(pth1+root)
data_resufled = read_milan_xyy_4col(data_names )


mossaic_keys = [['upper left', 'upper right'],
                ['lower left', 'lower right']]
fig, axes = plt.subplot_mosaic(mossaic_keys, constrained_layout=True,
                               gridspec_kw={'hspace': 0, 'wspace': 0})
fig.subplots_adjust(hspace=0)
read_plot_milan_ross_lags(data_resufled, axes)
prepare_figue_lag(fig, axes)
plt.savefig(pth1+"fig1_ross_lags.svg")

plt.show()

