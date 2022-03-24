
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


def plot_comoludogram(scales, phasPhas_TE, phasAmp_matrix):
    
    fig = plt.figure(figsize=(15, 7.5))

    a = np.linspace(0, 1, phasPhas_TE.shape[0])
    b = np.linspace(0, 1, phasPhas_TE.shape[1])

    xaxe, yaxe = np.meshgrid(a, b)

    gs = gridspec.GridSpec(1, 2)
    gs.update(left=0.05, right=0.95, hspace=0.3,
                top=0.95, bottom=0.05, wspace=0.15)
    axs = [gs[0, 0], gs[0, 1]]
    toplot = [phasPhas_TE.T, phasAmp_matrix.T]
    #toplot = [phase_TE, phase_amp_TE]
    tits = ['PHASE-PHASE CAUSALITY', 'PHASE-AMP CAUSALITY']
    labs = ['PHASE', 'AMP']
    max_mat = max(map(max, phasAmp_matrix))
    for ax, cont, tit, lab in zip(axs, toplot, tits, labs):
        ax = plt.subplot(ax)
        cs = ax.contourf(xaxe, yaxe, cont, levels=np.arange(
            0.0, max_mat, 0.00125), cmap=plt.cm.get_cmap("jet"), extend='max')
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
            plt.colorbar(cs)
            ax.grid()
            ax.set_ylabel("PERIOD %s [years]" % lab, size=23)
    #plt.savefig('plots/nino34-CESMhigh.eps', bbox_inches="tight")


for n in names:

    #scales, TransferEntropy = np.genfromtxt(n, delimiter=' ',
    #                                    dtype="f8,f8", unpack=True, usecols=[0, 1])
    TransferEntropy = np.genfromtxt(n, delimiter=' ',
                                            dtype="f8", unpack=True, usecols=[0])

    plot_TE_rows(n[-15:-4], TransferEntropy)

#fig = plt.figure('name')
#ax = plt.subplot(1, 1, 1) 
#names = read_texts(sys.argv[1]) 
#ax = plt.subplot(1, 1, 1) 
#plot_subfigures(ax, names)
# sys.argv[1]


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