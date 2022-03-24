
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
import pandas as pd
#import julia
#from julia.api import Julia
#jl = Julia(compiled_modules=False)


name = '../../package_CMI_prague/data/exp_raw/binfiles/Rossler_bin_0.000.bin'


#df = pd.read_csv('EURUSD.csv',sep='\t', index_col='Date')
#df = rsd.read_bin_bin_dataframe(name)
#df.sort_index(inplace=True)
#df = df.resample('W').last()
#sig =  np.array(df['x'][0:1000])


#131072 samples of the first components of the coupled Roessler systems, for 100 different epsilons from 0 to almost 0.25.
def read_raw_data(namem, start_index, end_index):
    usecols = [0,1]
    x, y= np.loadtxt(name,  comments='#', usecols = usecols, unpack = True, skiprows = start_index, max_rows= end_index )##, 

    return x, y


'''
def load_n_coupled_systems(n):
    ''''''loads the dat file n times and taques n series of corralations''''''
    coupling = np.linspace(0, max_coupling, 10)
    x_arr =[]
    y_arr =[]


    for c in coupling:

        start_index = int(c*100/max_coupling)*n_points
        end_index = start_index + n_points
        print (c, start_index , end_index  )
        x, y = read_raw_data(name, start_index, end_index)
        x_arr.append(x)
        y_arr.append(y)

    return x_arr, y_arr
'''
def write_coupled_bin(x, y, c):
    content = np.array([x, y])
    newFile = open('./binfiles/Rossler_bin_'+'{0:.3f}'.format(c)+'.bin', "wb")
    newFileByteArray = bytes(content)
    newFile.write(newFileByteArray)


def load_n_coupled_systems_read_once(num_couplings):
    '''loads the dat file and taques n series of corralations'''

    x, y = read_raw_data(name, 0, n_points*num_couplings)#num_couplings

    print ('total len', len(x))
    
    
    coupling = np.linspace(0, max_coupling, num_couplings-1)
    x_arr =[]
    y_arr =[]

    #7995392 8126464 131072

    for c in coupling:
        start_index = int(c*num_couplings/(max_coupling))*n_points
        end_index = start_index + n_points
        
        print ('number', int(c*num_couplings/(max_coupling)), 'coupling', c, start_index , end_index , len(x[start_index: end_index]) )
        #print (x[start_index], y[start_index])
        write_coupled_bin(x[start_index: end_index], y[start_index: end_index], c)
        
        x_arr.append(x[start_index: end_index])
        y_arr.append(y[start_index: end_index])

    return x_arr, y_arr

def read_bin_bin_dataframe(name):

    dt = np.dtype( '<f8')
    
    with open(name, 'rb') as f:
        b = f.read()
    np_data = np.frombuffer(b, dt)
    len_data = len(np_data)
    x = np_data[0:int(len_data/2)-1]
    y = np_data[int(len_data/2)+1:len_data]
    df = pd.DataFrame({'x':x, 'y':y})
    print (df)
    return df

n_points = 131072#000 #13107200
#name = '/Users/andreu/codes_TE/transfer_entropy_rossler_data/arosf11n00eps100raw.dat'
max_coupling = 0.25
#content = load_n_coupled_systems_read_once(100)
#name = '../data/exp_raw/binfiles/Rossler_bin_0.000.bin'
#read_bin_bin_dataframe(name)


#print ('len   '  , len(x_arr), len(y_arr),  len(x_arr[1]) )


#j = julia.Julia()
#a = j.include("/Users/andreu/codes_TE/Entropies_tests/src/Entropies_tests.jl")
#plt.plot(x[:end_index-1 ])
#plt.plot(y[:end_index-1 ])
#plt.show()