import pywt
import numpy as np
from numpy import random as rn

sampling_dt = 0.083
scales = np.array([1,2,3])
sig = rn.uniform(-1, 1, 555)

coeffs, freq_pywt = pywt.cwt(
                    sig,
                    scales=scales,
                    wavelet='cmor1.5-1.0',
                    sampling_period=sampling_dt,
                )
print(1/freq_pywt)