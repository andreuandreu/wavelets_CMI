from numpy import random as rd
import random as rd
import numpy as np



def circular_surrogate(data, min_shift= 1):
    '''Given a data series return it shifted by a random set of positions
     minimum and maximum shift shall be assured bonded by correlation length
     one doe not want to shift too little, or too much that it just goes around '''

    shift = int(rd.random()*(len(data)-2*min_shift)) + min_shift

    return np.roll(data, shift)
