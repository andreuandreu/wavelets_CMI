from numpy import random as rd
import random as rd
import numpy as np
#from pyunicorn.timeseries import surrogates as su



def circular_surrogate(data, min_shift= 1):
    '''Given a data series return it shifted by a random set of positions
     minimum and maximum shift shall be assured bonded by correlation length
     one doe not want to shift too little, or too much that it just goes around '''

    shift = int(rd.random()*(len(data)-2*min_shift)) + min_shift

    return np.roll(data, shift), shift



def many_surrogates(name, data, root = 'surr_', min_shift=1, n_surrogates = 111, txt = False):
    '''creates and stores a list of n surrogates given a timeseries '''
    #folder = './data/surrogates'

    #surrogates = np.zeros((n_surrogates, len(data)+1))
    #surrogates = np.array([])
    #surrogates = np.zeros(n_surrogates)
    surrogates = np.empty(n_surrogates, dtype=None)
    surrogates = []
    for i in range(n_surrogates):
        #if i % 11 == 1: print('surogate num', i)
        data_shifted, shift = circular_surrogate(data, min_shift)

        #surrogates[i] = np.append(data_shifted, shift)
        surrogates.append(np.append(data_shifted, shift))
        #np.append(surrogates, [data_shifted], axis=0)
        if txt:
            name_file = name + '_' + "{0:0>3}".format(i) + '.txt'
            np.savetxt(name ,    np.append(data_shifted, shift),
             
             
             
                     fmt='%.5e', newline='\n')
    
    name_file = name  + '.npy'
    np.save(name_file, surrogates)

data = np.array([np.arange(1,1111,4), np.arange(1,1111,4)])
#sur = su.Surrogates(data.astype('FLOATTYPE_t'))
#twin_sur = sur.twin_surrogates(data, dimension = 1, delay = 33, threshold = 4,
#                        min_dist=7)  
#print('im doing someting >))))', su)