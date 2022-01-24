from sklearn.feature_selection import mutual_info_classif
from entropy_estimators import continuous as eec
from sklearn.metrics import mutual_info_score
import numpy as np



def compure_MI_delays(dataX, dataY, delay):

    '''
    Given two data series of equal length, shift one of them up to a delay and compute their MI
    for two different methods
    '''

    MI_entEst_mnv = []
    MI_entEst = []
    MI_entEst_naive = []


    
    for d in range(1, delay):
        #mutual_info_classif(d1, d2)
        MI = eec.get_mi_mvn(dataX[:-d], dataY[d:]) #multivariate gausian distribution
        MI_entEst_mnv.append(MI)
        print(MI)
        MI_entEst.append(eec.get_mi(dataX[:-d], dataY[d:], k=16))
        MI_entEst_naive.append(eec.get_mi(
            dataX[:-d], dataY[d:], k=16, estimator='naive'))


    return np.array(MI_entEst), np.array(MI_entEst_mnv),  np.array(MI_entEst_naive)
