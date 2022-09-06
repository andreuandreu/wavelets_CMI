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

def filter_data( Z, cut = 2.0, filter = 3,  percentage = 0.66):

    
    #thresholdZ = max(Z) - (max(Z)-min(Z)) * percentage
    filtered_z = scipy.ndimage.zoom(Z, filter)#griddata(Z,  method='nearest')#gaussian_filter(Z, sigma=filter)#
    
    masked_array = np.ma.masked_where( (filtered_z > -cut) & (filtered_z < cut) , filtered_z)
    
    return masked_array

def Z_scoring(x):

    Zscore =  (np.array(x) - np.mean(x)) / np.std(x)
    return Zscore

def prepare_figure_comoludogram(fig, axes):

    left = 0.01
    width = 0.7
    bottom  = 0.01
    height = 0.8
    right = left + width
    top = bottom + height

    fig.text(0.86, 0.9, 'Z-score', va='center', fontsize=11)

    fig.text(0.45, 0.04, 'Pha[Hz]', va='center', fontsize=11)
    fig.text(0.048, 0.5, 'Amp[Hz]', va='center', rotation='vertical', fontsize=11)   

    fig.text(0.25, 0.9, r'Pha $\rightarrow$ Amp', va='center', fontsize=13)
    fig.text(0.65, 0.9, r'Amp $\rightarrow$ Pha', va='center', fontsize=13)

    fig.text(0.015, 0.2, 'PP-IC', va='center', rotation='vertical', fontsize=13)
    fig.text(0.015, 0.5, 'Im-IC', va='center', rotation='vertical', fontsize=13)
    fig.text(0.015, 0.76, 'Sch-IC', va='center', rotation='vertical', fontsize=13)


def prepare_figure_victor_como(fig, axes):


    fig.text(0.85, 0.92, 'Z-score', va='center', fontsize=11)
    fig.text(0.85,0.89, 'PAC', va='center', fontsize=11)
    fig.text(0.85, 0.09,  'APC', va='center', fontsize=11)

    fig.text(0.45, 0.04, 'Pha[Hz]', va='center', fontsize=11)
    #fig.text(0.06, 0.5, 'Amp[Hz]', va='center', rotation='vertical', fontsize=11)   

    #fig.text(0.02, 0.2, 'PP-IC', va='center', rotation='vertical', fontsize=13)
    #fig.text(0.02, 0.5, 'Im-IC', va='center', rotation='vertical', fontsize=13)
    #fig.text(0.02, 0.76, 'Sch-IC', va='center', rotation='vertical', fontsize=13)


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


def read_plot_all_milan_compludograms(data_names, fig, axes):

    #milan tags for identification data and direction
    #### v1, 2, 3 variables (Sch-IC, Im-IC, PP-IC) ###
    signame = { 'v1': 'Sch-IC',
                'v2': 'Im-IC',  
                'v3': 'PP-IC' }    
    direction =  { 'AAph': 'Amp - Pha',
                'phAA': 'Pha - Amp' }

    #place the plots in Milan's order, for now
    data_names = [data_names[3], data_names[5],
                 data_names[4],  data_names[1],
                 data_names[2], data_names[0]]
    maxZ = 1
    maxZs = []
    for key, n in zip( axes, data_names):
        dataX, dataY, dataZ = read_milan_xyz(n)
        
        if 'phAA' in n:
            maxZ =  np.max(dataZ)

        if 'AAph' in n:
            dataZ[0] = maxZ
            maxZs = np.append(maxZs, maxZ)
            divider = make_axes_locatable(axes[key])
            cax = divider.append_axes('right', size='5%', pad=0.05)
            fig.colorbar(im, cax=cax, orientation='vertical' )
           
        im = plot_milan_xyz(axes[key], dataX, dataY, dataZ)
        axes[key].set_transform(axes[key].transAxes)
    
    prepare_figure_comoludogram(fig, axes)
    plt.savefig(pth3+"fig3_milan_como_Z-score_smuth12.svg")
    
    return maxZs


def read_plot_victor_compludograms_by_Zvalue(name, fig, axes, maxZs):

    data = rmf.read_mat_victor(name)
    #data = [data, data]

    print('shape', np.shape(data), type (data))

    filter = 22
    for key, n, maxZ in zip( axes, data, maxZs):
        
        print('nnnn', n)
        dataX = data[n]['x']
        dataY = data[n]['y']
        dataZ = data[n]['PSI']

        Zscore  =  Z_scoring(dataZ)
        #,np.ma.masked_where(Zscore < 1.0, Zscore)
        #maxZ = np.max(dataZ)
        #dataZ[0] = maxZ

        #zoom and smooth the z matrix, reescale the x,y axes acordingly 
        lin = np.arange(len(dataX))
        Xinterp = scipy.interpolate.interp1d(lin, dataX)
        lin2 = np.linspace(0, len(dataX)-1 , num=len(dataX)*filter, endpoint=True ) 
        axeXpoints = Xinterp(lin2)

        lin = np.arange(len(dataY))
        Yinterp = scipy.interpolate.interp1d(lin, dataY)
        lin2 = np.linspace(0, len(dataY)-1 , num=len(dataY)*filter, endpoint=True ) 
        axeYpoints = Yinterp(lin2)

        #select the x interval of the plot
        cutXlow = np.where(axeXpoints > 5 )[0][0]
        cutXhigh = np.where(axeXpoints < 15 )[-1][-1]
        

        if 'left' in key:
            
            filtered_z = scipy.ndimage.zoom(Zscore, filter)
            masked_array = np.ma.masked_where( (filtered_z > -cut) , filtered_z)
            cutedZ = masked_array[:, cutXlow:cutXhigh]
            cutedZ[0, 0] = np.min(cutedZ )
            im = plot_victor_como(axes[key], axeXpoints[cutXlow:cutXhigh], axeYpoints, -cutedZ )
        
        if 'right' in key:
            filtered_z = scipy.ndimage.zoom(Zscore, filter)
            masked_array = np.ma.masked_where( (filtered_z < cut) , filtered_z)
            cutedZ = masked_array[:, cutXlow:cutXhigh]
            #cutedZ[0, 1] = np.max(cutedZ )
            im = plot_victor_como(axes[key], axeXpoints[cutXlow:cutXhigh], axeYpoints, cutedZ )



        

        divider = make_axes_locatable(axes[key])
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, orientation='vertical' )
           
        axes[key].set_transform(axes[key].transAxes)
        axes[key].set_yticks([])
    
    prepare_figure_victor_como(fig, axes)    
    plt.savefig(pth3+"fig3_victor_como_Z-score_smuth22.svg")


def read_plot_all_victor_compludograms(name, fig, axes, maxZs):

    data = rmf.read_mat_victor(name)

    filter = 22
    for key, n, maxZ in zip( axes, data, maxZs):
        
        dataX = data[n]['x']
        dataY = data[n]['y']
        dataZ = data[n]['PSI']

        Zscore  =  Z_scoring(dataZ)
        filteredZ = filter_data( Zscore, cut = 2.0, filter = filter)#,np.ma.masked_where(Zscore < 1.0, Zscore)
        #maxZ = np.max(dataZ)
        #dataZ[0] = maxZ

        #zoom and smooth the z matrix, reescale the x,y axes acordingly 
        lin = np.arange(len(dataX))
        Xinterp = scipy.interpolate.interp1d(lin, dataX)
        lin2 = np.linspace(0, len(dataX)-1 , num=len(dataX)*filter, endpoint=True ) 
        axeXpoints = Xinterp(lin2)

        lin = np.arange(len(dataY))
        Yinterp = scipy.interpolate.interp1d(lin, dataY)
        lin2 = np.linspace(0, len(dataY)-1 , num=len(dataY)*filter, endpoint=True ) 
        axeYpoints = Yinterp(lin2)

        #select the x interval of the plot
        cutXlow = np.where(axeXpoints > 5 )[0][0]
        cutXhigh = np.where(axeXpoints < 15 )[-1][-1]
        cutedZ = filteredZ[:, cutXlow:cutXhigh]

        #keep the color scale simmetric
        cutedZ[0, 0] = -np.min(cutedZ )
        cutedZ[0, 1] = -np.max(cutedZ )

        im = plot_victor_como(axes[key], axeXpoints[cutXlow:cutXhigh], axeYpoints, cutedZ )

        divider = make_axes_locatable(axes[key])
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, orientation='vertical' )
           
        axes[key].set_transform(axes[key].transAxes)
        axes[key].set_yticks([])
    
    prepare_figure_victor_como(fig, axes)    
    plt.savefig(pth3+"fig3_victor_como_Z-score_smuth22.svg")

def plot_victor_como(ax, dataX, dataY, dataZ):

   
    cmap = plt.get_cmap('PRGn')
    cmap.set_bad(color='white')
   
    im = ax.contourf( dataX, dataY,  dataZ,
                      extend='both', cmap=cmap)#contourlevels, 

    return(im)
        
            

def read_milan_xyz(name):

    data = np.genfromtxt(name, dtype=("f8", "f8" , "f8"), unpack=True, usecols=[0, 1, 3])

    dataX = data[0]
    dataY = data[1]
    dataZ = data[2]

    return dataX, dataY, dataZ

    

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




cut = 2.0 # cut for the significance of the Z score, dubious but mila's choice. 

'''FIGURE Comoludogram Milan'''
pth3 ="../../Desktop/Dropbox/transfer_inormation_prague/plots/fig3_comoludogram/"
root = "*.xyz"
data_names = glob(pth3+root)

mossaic_keys = [['upper left', 'upper right'],
                ['middle left', 'middle right'],
                ['lower left', 'lower right']]

hight = 7.8
fig, axes = plt.subplot_mosaic(mossaic_keys, sharex=True, sharey=True,
                              figsize=(6.1, hight), gridspec_kw={'hspace': 0, 'wspace': 0})
maxZs = read_plot_all_milan_compludograms(data_names, fig, axes)


'''FIGURE Comoludogram Victor'''
filename_mat = "../../Desktop/Dropbox/transfer_inormation_prague/data/imput/CFD_rats/CFD_S10.mat"

fig, axes = plt.subplot_mosaic(mossaic_keys, sharex=True, sharey=True,
                              figsize=(6.1, hight), gridspec_kw={'hspace': 0, 'wspace': 0})
read_plot_victor_compludograms_by_Zvalue(filename_mat, fig, axes, maxZs)

mossaic_keys = [['upper'],
                ['middle'],
                ['lower']]
fig, axes = plt.subplot_mosaic(mossaic_keys, sharex=True, sharey=True,
                              figsize=(3.1, hight), gridspec_kw={'hspace': 0, 'wspace': 0})
read_plot_all_victor_compludograms(filename_mat, fig, axes, maxZs)





plt.show()

