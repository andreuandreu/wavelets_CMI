import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
from glob import glob
from matplotlib import cm
import matplotlib as mpl


def sin_waves(amp = 1 ):

    x = np.arange(0, num_waves*2*np.pi)
    y = amp*np.sin(x)

    return x, y


x, y = sin_waves()

fig, ax = plt.subplots()

ax.plot(x, y)

plt.show()