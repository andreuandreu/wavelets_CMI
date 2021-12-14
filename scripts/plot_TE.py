
import pandas as pd
import matplotlib.pyplot as plt
import sys
import numpy as np
from scipy import signal
import Mutual_information_tools as mit

#python ./scripts/plot_TE.py ./data/output/corr_TE_rossler_sim_GF-NF_Hz-cmor1_amp-amp/Prob-est_VisFreq_b_bin-150_eDim-3_MI_each-tau.txt
def read_TransferEntropy_julia(name):
    ''' read formated file'''
    scales, TE = np.genfromtxt(name ,delimiter=' ' ,
                               dtype="f8,f8", unpack=True, usecols=[0, 1])

    return scales, TE

def plot_TE_rows(name,  TransferEntropy):

    fig = plt.figure('name')


    ax = plt.subplot(1, 1, 1)
    t = (TransferEntropy -np.mean(TransferEntropy) ) / np.std(TransferEntropy)
    ax.plot(t, label=name)
    #ax.set_xscale('log')
    ax.legend(frameon=False, fontsize = 10 , loc = 'lower left')

names = sys.argv[:]

for n in names:
    print('aaaaaaa', n)

    #scales, TransferEntropy = np.genfromtxt(n, delimiter=' ',
    #                                    dtype="f8,f8", unpack=True, usecols=[0, 1])
    TransferEntropy = np.genfromtxt(n, delimiter=' ',
                                            dtype="f8", unpack=True, usecols=[0])

    plot_TE_rows(n[-15:-4], TransferEntropy)

name_dataX = './data/output/rossler_phase_Nska_71Hz_niko_cmor1.5-1.0_pha.npy'
name_dataY = './data/output/rossler_phase_Nska_71Hz_niko_cmor1.5-1.0_pha.npy'
dataX = np.load(name_dataX )
dataY = np.load(name_dataY)


print('shape, size data', dataY.shape, dataX.size)
print('shape, size data', dataX[0].shape, dataX[0].size)
MI_entEst, MI_entEst_mnv, MI_score, MI_entEst_naive = mit.compure_MI_delays(dataX[0], dataY[0], 35)

                                    #    discrete_features=True) ) # mutual information of 0.69, expressed in nats

MI_score = (MI_score - np.mean(MI_score)) / \
    np.std(MI_score)
MI_entEst = (MI_entEst - np.mean(MI_entEst)) / \
    np.std(MI_entEst)
MI_entEst_mnv = (MI_entEst_mnv - np.mean(MI_entEst_mnv)) / \
    np.std(MI_entEst_mnv)
MI_entEst_naive = (MI_entEst_naive - np.mean(MI_entEst_naive)) / \
    np.std(MI_entEst_naive)

plt.plot(MI_entEst, label = 'Krasov est k=16')
#plt.plot(MI_entEst_mnv, label= 'multivariate Gaussian distribution')
plt.plot(MI_entEst_naive, label='Kozachenko-Leonenko k=16')


plt.show()