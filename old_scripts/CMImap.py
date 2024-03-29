"""
Computes MI/CMI analysis of NINO34 data.
This is the main script
"""


import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from datetime import date, datetime
from dateutil.relativedelta import relativedelta
from matplotlib.ticker import MultipleLocator, FuncFormatter
import cPickle
import sys
import matplotlib.gridspec as gridspec
import scipy.io as sio

COMPUTE = True # if True, the map will be evaluated, if False, it will be drawn
CMIP5model = None # None for data or name of the model + _ + number of TS as more time series is available
use_PRO_model = True

def read_rossler(fname):
    import csv
    r = {}
    f = open(fname, 'r')
    reader = csv.reader(f, lineterminator = "\n")
    eps = None
    for row in reader:
        if 'eps1' in row[0]:
            if eps is not None:
                r[eps] = np.array(r[eps])
            eps = float(row[0][8:])
            r[eps] = []
        if not "#" in row[0]:
            n = row[0].split()
            r[eps].append([float(n[0]), float(n[1])])
    r[eps] = np.array(r[eps])
    f.close()

    return r


if COMPUTE:
    import pyclits.wavelet_analysis as wvlt
    import pyclits.mutual_inf as MI
    from pyclits.geofield import DataField
    from pyclits.data_loaders import load_enso_index
    from pyclits.surrogates import SurrogateField
    from multiprocessing import Process, Queue
    # import cmi_estimation as cme


def load_enso_SSTs(num_ts = None, PROmodel = False, EMRmodel = None, DDEmodel = None):
    # load enso SSTs
    print("[%s] Loading monthly ENSO SSTs..." % str(datetime.now()))
    # enso_raw = np.loadtxt("nino34m13.txt") # length x 2 as 1st column is continuous year, second is SST in degC
    # enso = DataField()

    # enso = load_enso_index("nino34raw.txt", '3.4', date(1900, 1, 1), date(2001, 1, 1))
    enso = DataField()
    enso.data = np.loadtxt("nino34_CESMhigh.dat")
    enso.create_time_array(date_from = date(1900, 1, 1), sampling = 'm')
    # fname = "conceptualRossler1:2monthlysampling_100eps0-0.25.dat"
    # r = read_rossler(fname)
    # x, y = r[0.0707][20000:52768, 0], r[0.0707][20000:52768, 1] # 
    # x, y = r[0.1465][20000:21024, 0], r[0.1465][20000:21024, 1] # x is annual, y is biennal
    # print("rossler data read")

    # add noise
    # x += np.random.normal(0, scale = 0.1*np.std(x, ddof = 1), size = (x.shape[0],))
    # y += np.random.normal(0, scale = 0.1*np.std(y, ddof = 1), size = (y.shape[0],))
    # print("added noise to rossler data...")
    # exa = np.loadtxt("ExA-comb-mode-20CR-1900-2010-PC2-stand.txt")
    # enso.data = exa.copy()
    # enso.data = enso_raw[:, 1]
    # enso.data = np.zeros((1200,))
    # if '4k' in EMRmodel:
    #     enso.data = np.zeros((4096,))
    # elif '8k' in EMRmodel:  
    #     enso.data = np.zeros((8192,))
    # elif '16k' in EMRmodel:
    #     enso.data = np.zeros((16384,))
    # elif '32k' in EMRmodel:
    #     enso.data = np.zeros((32768,))
    # enso.data = x.copy()
    # enso.create_time_array(date(1900, 1, 1), sampling = 'm')
    print enso.data.shape
    print enso.get_date_from_ndx(0), enso.get_date_from_ndx(-1)

    # enso.select_date(date(1884,1,1), date(2014,1,1))

    # if CMIP5model is not None and num_ts is not None:
    #     fname = CMIP5model + '.txt'
    #     model = np.loadtxt('N34_CMIP5/' + fname)
    #     enso.data = model[:, num_ts]
    # with open("SST-EOFanalysis-inputPCs.bin", "rb") as f:
    #     pcs = cPickle.load(f)
    # pcs = pcs['PC']
    # if num_ts != 2:
    #     enso.data = pcs[0, :]
    # elif num_ts == 2:
    #     enso.data = pcs[3, :]

    if PROmodel:
        print("[%s] Integrating PRO model which will be used instead of ENSO SSTs..." % (str(datetime.now())))
        from parametric_recharge_oscillator import ENSO_PROmodel
        PROmodel_enso = ENSO_PROmodel(length = enso.data.shape[0], daily = False, damped = False, ENSOperiod = 3.75, modulation = 2., lambda0 = 0.4)
        PROmodel_enso.integrate_PROmodel()
        enso.data = PROmodel_enso.data.copy()

    if DDEmodel is not None:
        print("[%s] Integrating DDE model which will be used instead of ENSO SSTs..." % (str(datetime.now())))
        from enso_dde_model import enso_dde_model
        kappa, tau, b, subsample = DDEmodel
        t, h = enso_dde_model(kappa = kappa, tau = tau, b = b, nyears = enso.time.shape[0]//12+10, subsample_to_monthly = subsample)
        enso.data = h.copy()

    if EMRmodel is not None:
        print("[%s] Loading EMR (%s) simulated syntethic ENSO time series..." % (str(datetime.now()), EMRmodel))
        # raw = sio.loadmat("Nino34-%s.mat" % (EMRmodel))['N34s']
        raw = sio.loadmat("Nino34-linear-10PCsel-L3-seasonal-surrs.mat")['N34s']
        enso.data = raw[num_ts, -1332:] # same length as nino3.4 data
        # if '4k' in EMRmodel:
        #     enso.data = raw[-4096:, num_ts].copy()
        # elif '8k' in EMRmodel:  
        #     enso.data = raw[-8192:, num_ts].copy()
        # elif '16k' in EMRmodel:
        #     enso.data = raw[-16384:, num_ts].copy()
        # elif '32k' in EMRmodel:
        #     enso.data = raw[-32768:, num_ts].copy()

        # if num_ts%3 == 0:
        #     dat = enso.get_date_from_ndx(12423)
        # elif num_ts%3 == 1:
        #     dat = enso.get_date_from_ndx(2434)
        # elif num_ts%3 == 2:
        #     dat = enso.get_date_from_ndx(7354)
        
    # enso.get_data_of_precise_length(length = 1200, end_date = enso.get_date_from_ndx(-1), apply_to_data = True)


    if NUM_SURR > 0:
        a = list(enso.get_seasonality(detrend = False))
        enso_sg = SurrogateField()

        # _, _, idx = enso.get_data_of_precise_length(length = 1024, end_date = date(2014, 1, 1), COPY = False)
        enso_sg.copy_field(enso)
        # enso_sg.data = enso_sg.data[idx[0]:idx[1]]

        enso.return_seasonality(a[0], a[1], None)

        # a[0] = a[0][idx[0]:idx[1]]
        # a[1] = a[1][idx[0]:idx[1]]

    # select 1024 data points
    # enso.get_data_of_precise_length(length = 1024, end_date = date(2014, 1, 1), COPY = True)
    print("[%s] Data loaded with shape %s" % (str(datetime.now()), enso.data.shape))

    return enso, enso_sg, a#, y

def phase_diff(ph1, ph2):
    ph = ph1 - ph2
    ph[ph < -np.pi] += 2*np.pi

    return ph


WVLT_SPAN = [5, 96]  # unit is month, 96
NUM_SURR = 100
WRKRS = 20
# BINS = 4
bins_list = [4]

# CMIP5model = 'N34_CanESM2_0'# None for data or name of the model + _ + number of TS as more time series is available
# CMIP5models = ['N34_CanESM2', 'N34_GFDLCM3', 'N34_GISSE2Hp1', 'N34_GISSE2Hp2', 'N34_GISSE2Hp3', 'N34_GISSE2Rp1']
# CMIP5models += ['N34_GISSE2Rp2', 'N34_GISSE2Rp3', 'N34_HadGem2ES', 'N34_IPSL_CM5A_LR', 'N34_MIROC5', 'N34_MRICGCM3']
# CMIP5models += ['N34_CCSM4', 'N34_CNRMCM5', 'N34_CSIROmk360']
#CMIP5models = [[50., 0.42, 1., True], [5., 0.65, 1., True], [50., 0.508, 1., True], [11., 0.56, 1.4, True]]
CMIP5models = ['delay-quad-10PCsel-L3-seasonal-d:7mon-k:50']

if COMPUTE:
    for BINS in bins_list:
        print("[%s] Computing using %d bins" % (str(datetime.now()), BINS))
        for CMIP5model in CMIP5models:
            # fname = CMIP5model + '.txt'
            # model = np.loadtxt('N34_CMIP5/' + fname)
            # model_count = model.shape[1]
            model_count = 1
            # exa = np.loadtxt("ExA-comb-mode-20CR-1900-2010-PC2-stand.txt")[-1332:]
            # exa = np.loadtxt("PC1-wind-comb-mode-20CR-1900-2010-stand.txt")[-1332:]
            # CMIP5model = None

            for num_ts in range(model_count):
            # for num_ts in model_count:

                # print("[%s] Evaluating %d. time series of %s model data... (%d out of %d models)" % (str(datetime.now()), 
                    # num_ts, CMIP5model, CMIP5models.index(CMIP5model)+1, len(CMIP5models)))

                enso, enso_sg, seasonality = load_enso_SSTs(num_ts, PROmodel = False, EMRmodel = None, DDEmodel = None)
                #exa = sio.loadmat('Nino34-SST-x-wind-40PCsel.mat')['wind_sim'][:, 1, num_ts]
                # if num_ts == 0:
                #     exa = enso.data.copy() # 1 -> 1
                # elif num_ts == 1:
                #     exa = pcs[3, :] # 1 -> 4
                # elif num_ts == 2:
                #     exa = pcs[0, :] # 4 -> 1
                # qbo = np.loadtxt("QBOindex-jan1948-dec2015.txt")
                ## DATA
                #prepare result matrices
                # 1st PCA of spatial QBO
                # qbo = np.loadtxt("/Users/nikola/work-ui/data/QBOspatial-2-PCA-jan1956-dec2015.txt")[:, 0]
                scales = np.arange(WVLT_SPAN[0], WVLT_SPAN[-1] + 1, 1)
                phase_phase_coherence = np.zeros((scales.shape[0], scales.shape[0]))
                # phase_phase_coherence = np.zeros((scales.shape[0], 1))
                phase_phase_CMI = np.zeros_like(phase_phase_coherence)
                phase_amp_MI = np.zeros_like(phase_phase_coherence)
                phase_amp_condMI = np.zeros_like(phase_phase_coherence)
                phase_phase_coherence_knn = np.zeros((scales.shape[0], scales.shape[0]))
                # phase_phase_coherence_knn = np.zeros((scales.shape[0], 1))
                phase_phase_CMI_knn = np.zeros_like(phase_phase_coherence_knn)
                phase_amp_MI_knn = np.zeros_like(phase_phase_coherence_knn)
                phase_amp_condMI_knn = np.zeros_like(phase_phase_coherence_knn)

                enso.center_data()
                # qbo -= np.nanmean(qbo, axis = 0)
                # y_ts -= np.nanmean(y_ts, axis = 0)
                print enso.data.shape#, y_ts.shape

                for i in range(phase_phase_coherence.shape[0]):
                    enso.wavelet(period = scales[i], period_unit = 'm', cut = 1)
                    phase_i = enso.phase.copy()
                    
                    for j in range(phase_phase_coherence.shape[1]):
                        enso.wavelet(period = scales[j], period_unit = 'm', cut = 1)
                        phase_j = enso.phase.copy()
                        
                        # ###################
                        # # continuous phase
                        # for t in range(phase_j.shape[0] - 1):
                        #     if np.abs(phase_j[t+1] - phase_j[t]) > 1.:
                        #         phase_j[t+1:] += 2*np.pi
                        
                        # # phase fluctuations
                        # omega = 2*np.pi / (12 * 1) # annual
                        # ph0 = phase_j[0]
                        # for t in range(phase_j.shape[0]):
                        #     phase_j[t] -= (ph0 + omega*t)
                        # ###################

                        amp_j = enso.amplitude.copy()

                        phase_phase_coherence[i, j] = MI.mutual_information(phase_i, phase_j, algorithm = 'EQQ2', bins = BINS)
                        phase_phase_coherence_knn[i, j] = MI.knn_mutual_information(phase_i, phase_j, k = 64, dualtree = True)

                        # conditional mutual inf -- avg over lags 1 - 6 months
                        CMI = []
                        CMI_knn = []
                        for tau in range(1, 7):
                            CMI.append(MI.cond_mutual_information(phase_i[:-tau], phase_diff(phase_j[tau:], phase_j[:-tau]), 
                                phase_j[:-tau], algorithm = 'EQQ2', bins = BINS))
                            CMI_knn.append(MI.knn_cond_mutual_information(phase_i[:-tau], phase_diff(phase_j[tau:], phase_j[:-tau]), 
                                phase_j[:-tau], k = 64, dualtree = True))
                        phase_phase_CMI[i, j] = np.mean(np.array(CMI))
                        phase_phase_CMI_knn[i, j] = np.mean(np.array(CMI_knn))

                        phase_amp_MI[i, j] = MI.mutual_information(phase_i, amp_j, algorithm = 'EQQ2', bins = BINS)
                        phase_amp_MI_knn[i, j] = MI.knn_mutual_information(phase_i, amp_j, k = 64, dualtree = True)

                        CMI2 = []
                        CMI2_knn = []
                        eta = np.int(scales[i] / 4)
                        for tau in range(1, 7): # possible 1-31
                            x, yts, z = MI.get_time_series_condition([phase_i, np.power(amp_j,2)], tau = tau, dim_of_condition = 3, eta = eta)
                            CMI2.append(MI.cond_mutual_information(x, yts, z, algorithm = 'GCM', bins = BINS))
                            CMI2_knn.append(MI.knn_cond_mutual_information(x, yts, z, k = 64, dualtree = True))
                        phase_amp_condMI[i, j] = np.mean(np.array(CMI2))
                        phase_amp_condMI_knn[i, j] = np.mean(np.array(CMI2_knn))

                print("[%s] Analysis on data done." % str(datetime.now()))


                def _coh_cmi_surrs(sg, a, sc, jobq, resq):
                    mean, var, _ = a
                    # while s is not None:
                    while True:
                        s = jobq.get()
                        if s is None:
                            break
                        sg.construct_fourier_surrogates(algorithm = 'FT')
                        sg.add_seasonality(mean, var, None)

                        # sg.surr_data = s.copy()

                        coh = np.zeros((sc.shape[0], sc.shape[0]))
                        # coh = np.zeros((sc.shape[0], 1))
                        cmi = np.zeros_like(phase_phase_coherence)
                        ph_amp_MI = np.zeros_like(phase_phase_coherence)
                        ph_amp_CMI = np.zeros_like(phase_phase_coherence)
                        coh_knn = np.zeros((sc.shape[0], sc.shape[0]))
                        # coh_knn = np.zeros((sc.shape[0], 1))
                        cmi_knn = np.zeros_like(phase_phase_coherence_knn)
                        ph_amp_MI_knn = np.zeros_like(phase_phase_coherence_knn)
                        ph_amp_CMI_knn = np.zeros_like(phase_phase_coherence_knn)

                        sg.center_surr()

                        for i in range(coh.shape[0]):
                            sg.wavelet(period = sc[i], period_unit = 'm', cut = 1)
                            phase_i = sg.phase.copy()
                            
                            for j in range(coh.shape[1]):
                                sg.wavelet(period = sc[j], period_unit = 'm', cut = 1)
                                phase_j = sg.phase.copy()

                                # ###################
                                # # continuous phase
                                # for t in range(phase_j.shape[0] - 1):
                                #     if np.abs(phase_j[t+1] - phase_j[t]) > 1.:
                                #         phase_j[t+1:] += 2*np.pi
                                
                                # # phase fluctuations
                                # omega = 2*np.pi / (12 * 1) # annual
                                # ph0 = phase_j[0]
                                # for t in range(phase_j.shape[0]):
                                #     phase_j[t] -= (ph0 + omega*t)
                                # ###################

                                amp_j = sg.amplitude.copy()

                                coh[i, j] = MI.mutual_information(phase_i, phase_j, algorithm = 'EQQ2', bins = BINS)
                                coh_knn[i, j] = MI.knn_mutual_information(phase_i, phase_j, k = 64, dualtree = True)

                                # conditional mutual inf -- avg over lags 1 - 6 months
                                CMI_temp = []
                                CMI_temp_knn = []
                                for tau in range(1, 7):
                                    CMI_temp.append(MI.cond_mutual_information(phase_i[:-tau], phase_diff(phase_j[tau:], phase_j[:-tau]), 
                                        phase_j[:-tau], algorithm = 'EQQ2', bins = BINS))
                                    CMI_temp_knn.append(MI.knn_cond_mutual_information(phase_i[:-tau], phase_diff(phase_j[tau:], phase_j[:-tau]), 
                                        phase_j[:-tau], k = 64, dualtree = True))
                                cmi[i, j] = np.mean(np.array(CMI_temp))
                                cmi_knn[i, j] = np.mean(np.array(CMI_temp_knn))

                                ph_amp_MI[i, j] = MI.mutual_information(phase_i, amp_j, algorithm = 'EQQ2', bins = BINS)
                                ph_amp_MI_knn[i, j] = MI.knn_mutual_information(phase_i, amp_j, k = 64, dualtree = True)

                                CMI2 = []
                                CMI2_knn = []
                                eta = np.int(scales[i] / 4)
                                for tau in range(1, 7): # possible 1-31
                                    x, yts, z = MI.get_time_series_condition([phase_i, np.power(amp_j,2)], tau = tau, dim_of_condition = 3, eta = eta)
                                    CMI2.append(MI.cond_mutual_information(x, yts, z, algorithm = 'GCM', bins = BINS))
                                    CMI2_knn.append(MI.knn_cond_mutual_information(x, yts, z, k = 64, dualtree = True))
                                ph_amp_CMI[i, j] = np.mean(np.array(CMI2))
                                ph_amp_CMI_knn[i, j] = np.mean(np.array(CMI2_knn))

                        resq.put((coh, cmi, ph_amp_MI, ph_amp_CMI, coh_knn, cmi_knn, ph_amp_MI_knn, ph_amp_CMI_knn))
                        # resq.put((coh, cmi, ph_amp_MI, ph_amp_CMI))


                ## SURROGATES
                if NUM_SURR > 0:
                    print("[%s] Analysing %d FT surrogates using %d workers..." % (str(datetime.now()), NUM_SURR, WRKRS))

                    # surrs = sio.loadmat("DimaKon-Nino34-ERM-linear-SSTA.mat")['sstn']
                    # surrs = sio.loadmat("10m-wind-20PCs-L3-model-surrs.mat")['ExA_mode']
                    # surrs = surrs[-1024:, :].copy()
                    
                    surr_completed = 0
                    surrCoherence = np.zeros(([NUM_SURR] + list(phase_phase_coherence.shape)))
                    surrCMI = np.zeros_like(surrCoherence)
                    surrPhaseAmp = np.zeros_like(surrCoherence)
                    surrPhaseAmpCMI = np.zeros_like(surrCoherence)
                    surrCoherence_knn = np.zeros(([NUM_SURR] + list(phase_phase_coherence_knn.shape)))
                    surrCMI_knn = np.zeros_like(surrCoherence_knn)
                    surrPhaseAmp_knn = np.zeros_like(surrCoherence_knn)
                    surrPhaseAmpCMI_knn = np.zeros_like(surrCoherence_knn)
                    jobq = Queue()
                    resq = Queue()
                    for i in range(NUM_SURR):
                        # jobq.put(surrs[:, i])
                        jobq.put(1)
                    for i in range(WRKRS):
                        jobq.put(None)

                    wrkrs = [Process(target=_coh_cmi_surrs, args = (enso_sg, seasonality, scales, jobq, resq)) for i in range(WRKRS)]
                    for w in wrkrs:
                        w.start()

                    while surr_completed < NUM_SURR:
                        coh, cmi, phAmp, phAmpCMI, coh_knn, cmi_knn, phAmp_knn, phAmpCMI_knn = resq.get()
                        # coh, cmi, phAmp, phAmpCMI = resq.get()
                        surrCoherence[surr_completed, :, :] = coh
                        surrCMI[surr_completed, :, :] = cmi
                        surrPhaseAmp[surr_completed, :, :] = phAmp
                        surrPhaseAmpCMI[surr_completed, :, :] = phAmpCMI
                        surrCoherence_knn[surr_completed, :, :] = coh_knn
                        surrCMI_knn[surr_completed, :, :] = cmi_knn
                        surrPhaseAmp_knn[surr_completed, :, :] = phAmp_knn
                        surrPhaseAmpCMI_knn[surr_completed, :, :] = phAmpCMI_knn
                        surr_completed += 1

                        if surr_completed % 20 == 0:
                            print("..%d/%d surrogate done.." % (surr_completed, NUM_SURR))

                    for w in wrkrs:
                        w.join()

                    print("[%s] %d surrogates done. Saving..." % (str(datetime.now()), NUM_SURR))


                # fname = ("CMImap%dbins3Dcond_GaussCorr_%sts%d.bin" % (BINS, CMIP5model, num_ts))
                if use_PRO_model:
                    fname = ("PROneutral-CMImap%dbins3Dcond_GaussCorr.bin" % (BINS))
                # fname = ("DDEmodel-k%.1f-tau:%.3f-b:%.1f-against%dFT.bin" % (CMIP5model[0], CMIP5model[1], CMIP5model[2], NUM_SURR))
                # fname = ("kNN-Nino34-obs_CMImap4bins3Dcond-vs-Dima.bin")
                # fname = ("conceptualRossler-NOsynch-1:2-monthlyEQQ.bin")
                # fname = ("conceptualRossler-eps=0.1465-synch-1:2-monthlyEQQ-and-kNN-added-noise.bin")
                # fname = ("qbo-to-nino34-1948-2015-monthly-EQQonly.bin")
                # fname = ("nino34-1870-1943-monthly-EQQandKNN-500FT-1-30avgCMI.bin")
                # fname = ("nino34-1900-2010-monthly-EQQandKNN-annual-ph-flucts.bin")
                # fname = ("SST-PCs-type%d_CMImap4bins3Dcond-against-500FT.bin" % (num_ts))
                # fname = ("kNN-PROdamped-3.75per_CMImap4bins3Dcond%d-against-500FT.bin" % (num_ts))
                # fname = ("PC1-wind-vs-ExA-comb-mode-as-x-vs-y_CMImap4bins3Dcond-500FT.bin")
                fname = ("nino34-CESMhigh.bin")
                with open(fname, 'wb') as f:
                    cPickle.dump({'phase x phase data' : phase_phase_coherence, 'phase CMI data' : phase_phase_CMI, 
                        'phase x phase surrs' : surrCoherence, 'phase CMI surrs' : surrCMI, 'phase x amp data' : phase_amp_MI,
                        'phase amp CMI data' : phase_amp_condMI, 'phase x amp surrs' : surrPhaseAmp, 'phase amp CMI surrs' : surrPhaseAmpCMI, 
                        'phase x phase data knn' : phase_phase_coherence_knn, 'phase CMI data knn' : phase_phase_CMI_knn, 
                        'phase x phase surrs knn' : surrCoherence_knn, 'phase CMI surrs knn' : surrCMI_knn, 'phase x amp data knn' : phase_amp_MI_knn,
                        'phase amp CMI data knn' : phase_amp_condMI_knn, 'phase x amp surrs knn' : surrPhaseAmp_knn, 'phase amp CMI surrs knn' : surrPhaseAmpCMI_knn}, 
                        f, protocol = cPickle.HIGHEST_PROTOCOL)
        print("[%s] All models done." % str(datetime.now()))



else:
    CMIP5models = ['linear-20PC-L3-seasonal']
    BINS = 4
    PUB = False
    # WVLT_SPAN = [5,93]
    for CMIP5model in CMIP5models:
        # fname = CMIP5model + '.txt'
        # model = np.loadtxt('N34_CMIP5/' + fname)
        # model_count = model.shape[1]
        if PUB:
            model_count = 1
        else:
            model_count = 1
        # CMIP5model = None
        scales = np.arange(WVLT_SPAN[0], WVLT_SPAN[-1] + 1, 1)
        overall_ph_ph = np.zeros((scales.shape[0], scales.shape[0]))
        overall_ph_ph_cmi = np.zeros_like(overall_ph_ph)
        overall_ph_amp = np.zeros_like(overall_ph_ph)
        overall_ph_amp_cmi = np.zeros_like(overall_ph_ph)

        # model_count = ['34']

        for num_ts in range(model_count):
            # fname = ("bins/python-model/Python-Nino34-%s_CMImap4bins3Dcond%d-against-basicERM.bin" % (CMIP5model, num_ts))
            # fname = ("bins/Nino%s-obs-vs-ExA-Sergey-reversed-comb-mode_CMImap4bins3Dcond-against-basicERM.bin" % (num_ts))
            # fname = ("bins/Sergey-Nino34-ERM-linear-SST-20PC-L3-multiplicative-5mon-snippets_CMImap4bins3Dcond%d-against-basicERM.bin" % (num_ts))
            fname = 'bins/nino34-CESMhigh.bin'
            CUT = slice(0,NUM_SURR)
            # version = 3
            if num_ts == 8:
                continue
            with open(fname, 'rb') as f:
                result = cPickle.load(f)

            phase_phase_coherence = result['phase x phase data knn']
            phase_phase_CMI = result['phase CMI data knn']
            surrCoherence = result['phase x phase surrs knn']
            surrCMI = result['phase CMI surrs knn']
            phase_amp_MI = result['phase x amp data knn']
            phase_amp_condMI = result['phase amp CMI data knn']
            surrPhaseAmp = result['phase x amp surrs knn']
            surrPhaseAmpCMI = result['phase amp CMI surrs knn']

            print num_ts, "loaded"

            res_phase_coh = np.zeros_like(phase_phase_coherence)
            res_phase_cmi = np.zeros_like(res_phase_coh)
            res_phase_amp = np.zeros_like(res_phase_coh)
            res_phase_amp_CMI = np.zeros_like(res_phase_coh)

            for i in range(res_phase_coh.shape[0]):
                for j in range(res_phase_coh.shape[1]):
                    res_phase_coh[i, j] = np.sum(np.greater(phase_phase_coherence[i, j], surrCoherence[:, i, j])) / np.float(surrCoherence.shape[0])
                    # res_phase_coh[i, j] = (phase_phase_coherence[i, j] - np.mean(surrCoherence[:, i, j])) / np.std(surrCoherence[:, i, j], ddof = 1)
                    res_phase_cmi[i, j] = np.sum(np.greater(phase_phase_CMI[i, j], surrCMI[:, i, j])) / np.float(surrCMI.shape[0])
                    res_phase_amp[i, j] = np.sum(np.greater(phase_amp_MI[i, j], surrPhaseAmp[:, i, j])) / np.float(surrPhaseAmp.shape[0])
                    res_phase_amp_CMI[i, j] = np.sum(np.greater(phase_amp_condMI[i, j], surrPhaseAmpCMI[:, i, j])) / np.float(surrPhaseAmpCMI.shape[0])

            overall_ph_ph += np.greater_equal(res_phase_coh, 0.95)
            overall_ph_ph_cmi += np.greater_equal(res_phase_cmi, 0.95)
            overall_ph_amp += np.greater_equal(res_phase_amp, 0.95)
            overall_ph_amp_cmi += np.greater_equal(res_phase_amp_CMI, 0.95)

            # a = np.random.rand(scales.shape[0], scales.shape[0]) + 0.5
            x, y = np.meshgrid(scales, scales)

            # fig, axs = plt.subplots(1, 2, figsize = (13,7))
            if PUB:
                fig = plt.figure(figsize=(15,7.5))
                gs = gridspec.GridSpec(1, 2)
                gs.update(left=0.05, right=0.95, hspace=0.3, top=0.95, bottom=0.05, wspace=0.15)
                axs = [gs[0,0], gs[0,1]]
                plot = [res_phase_cmi.T, res_phase_amp_CMI.T]
                # plot = [res_phase_coh.T, res_phase_cmi.T]
                tits = ['PHASE-PHASE CAUSALITY', 'PHASE-AMP CAUSALITY']
                # tits = ['PHASE SYNCHRONIZATION', 'PHASE-PHASE CAUSALITY']
                labs = ['PHASE', 'AMP']
                # labs = ['PHASE', 'PHASE']
                for ax, cont, tit, lab in zip(axs, plot, tits, labs):
                    ax = plt.subplot(ax)
                    cs = ax.contourf(x, y, cont, levels = np.arange(0.95, 1, 0.00125), cmap = plt.cm.get_cmap("jet"), extend = 'max')
                    # cs = ax.contourf(x, y, cont, levels = np.arange(4, 20, 0.125), cmap = plt.cm.get_cmap("jet"), extend = 'max')
                    ax.tick_params(axis='both', which='major', labelsize = 20)
                    ax.set_title(tit, size = 30)
                    ax.xaxis.set_major_locator(MultipleLocator(12))
                    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: int(x)/12))
                    ax.xaxis.set_minor_locator(MultipleLocator(6))
                    ax.yaxis.set_major_locator(MultipleLocator(12))
                    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: int(x)/12))
                    ax.yaxis.set_minor_locator(MultipleLocator(6))
                    ax.set_xlabel("PERIOD PHASE [years]", size = 23)
                    # plt.colorbar(cs)
                    ax.grid()
                    ax.set_ylabel("PERIOD %s [years]" % lab, size = 23)
                plt.savefig('plots/nino34-CESMhigh.eps', bbox_inches = "tight")
            else:
                fig = plt.figure(figsize=(15,15))
                gs = gridspec.GridSpec(2, 2)
                gs.update(left=0.05, right=0.95, hspace=0.3, top=0.95, bottom=0.05, wspace=0.15)
                i = 0
                axs = [gs[0,0], gs[0,1], gs[1,0], gs[1,1]]
                plot = [res_phase_coh.T, res_phase_cmi.T, res_phase_amp.T, res_phase_amp_CMI.T]
                tits = ['PHASE SYNCHRONIZATION', 'PHASE-PHASE CAUSALITY', 'PHASE x AMP MI', 'PHASE-AMP CAUSALITY']
                for ax, cont, tit in zip(axs, plot, tits):
                    ax = plt.subplot(ax)
                    cs = ax.contourf(x, y, cont, levels = np.arange(0.95, 1, 0.00125), cmap = plt.cm.get_cmap("jet"), extend = 'max')
                    ax.tick_params(axis='both', which='major', labelsize = 17)
                    ax.set_title(tit, size = 28)
                    ax.xaxis.set_major_locator(MultipleLocator(12))
                    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: int(x)/12))
                    ax.xaxis.set_minor_locator(MultipleLocator(6))
                    ax.yaxis.set_major_locator(MultipleLocator(12))
                    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: int(x)/12))
                    ax.yaxis.set_minor_locator(MultipleLocator(6))
                    ax.set_xlabel("annual period [years]", size = 20)
                    # plt.colorbar(cs)
                    ax.grid()
                    if i % 2 == 0:
                        ax.set_ylabel("biennal period [years]", size = 20)
                    else:
                        # fig.colorbar(cs, ax = ax, shrink = 0.5)
                        pass
                    i += 1

                # plt.savefig('plots/PNino34-%s-CMImap4bin%d-against-basicERM.png' % (CMIP5model, num_ts))
                # plt.savefig("plots/DDEmodel-k%.1f-tau:%.3f-b:%.1f-against%dFT.png" % (CMIP5model[0], CMIP5model[1], CMIP5model[2], NUM_SURR))
                # plt.savefig('plots/Nino%s-obs-vs-ExA-Sergey-reversed-comb-mode_CMImap4bins3Dcond-against-basicERM.png' % num_ts)
                # plt.savefig('plots/pro_knn/kNN-PROdamped-3.75per_CMImap4bins3Dcond%d-FT.png' % num_ts)
                plt.savefig('plots/nino34-CESMhigh.png')
            # plt.savefig('PROdamped-CMImap.png')
        # plt.savefig('test.png')

        if not PUB:
            fig = plt.figure(figsize=(15,15))
            gs = gridspec.GridSpec(2, 2)
            gs.update(left=0.05, right=0.95, hspace=0.3, top=0.95, bottom=0.05, wspace=0.15)
            i = 0
            axs = [gs[0,0], gs[0,1], gs[1,0], gs[1,1]]
            plot = [overall_ph_ph.T, overall_ph_ph_cmi.T, overall_ph_amp.T, overall_ph_amp_cmi.T]
            tits = ['PHASE SYNCHRONIZATION', 'PHASE-PHASE CAUSALITY', 'PHASE x AMP MI', 'PHASE-AMP CAUSALITY']
            labs = ['PHASE', 'PHASE', 'AMP', 'AMP']
            for ax, cont, tit, lab in zip(axs, plot, tits, labs):
                ax = plt.subplot(ax)
                cs = ax.contourf(x, y, cont, levels = np.arange(1, model_count+1, 1), cmap = plt.cm.get_cmap("jet"))
                # cs = ax.pcolormesh(x, y, cont, vmin = 1)
                ax.tick_params(axis='both', which='major', labelsize = 20)
                ax.set_title(tit, size = 30)
                ax.xaxis.set_major_locator(MultipleLocator(12))
                ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: int(x)/12))
                ax.xaxis.set_minor_locator(MultipleLocator(6))
                ax.yaxis.set_major_locator(MultipleLocator(12))
                ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: int(x)/12))
                ax.yaxis.set_minor_locator(MultipleLocator(6))
                ax.set_xlabel("PERIOD PHASE [years]", size = 23)
                plt.colorbar(cs)
                ax.grid()
                if i % 1 == 0:
                    ax.set_ylabel("PERIOD %s [years]" % lab, size = 23)
                else:
                    # fig.colorbar(cs, ax = ax, shrink = 0.5)
                    pass
                i += 1

            # plt.savefig('plots/PNino34-%s-CMImap4bin-overall.png' % (CMIP5model), bbox_inches = "tight")
            # plt.savefig('plots/pro_knn/kNN-PROdamped-3.75per_CMImap4bins3Dcond-overall.eps', bbox_inches = "tight")
            plt.savefig("ERMplot.eps", bbox_inches = "tight")
            # plt.savefig('plots/wind-x-sst-model/Nino-obs_CMImap4bins3Dcond-against-basicERM-overall.png', bbox_inches = "tight")

            print np.unique(overall_ph_ph)
