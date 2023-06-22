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


def read_milan_xyy_4col(names):

    data = np.empty((4, 2, 50))

    for i, n in enumerate(names):
        d = np.genfromtxt(
            n, dtype=("f8", "f8", "f8", "f8"), unpack=True, usecols=[1, 2, 3, 4]
        )
        data[2 * i] = d[0:2]
        data[2 * i + 1] = d[2:4]

    return data


def read_milan_xyy_2col(name):

    data = np.genfromtxt(name, dtype=("f8", "f8"), unpack=True, usecols=[1, 2])
    return data


def read_plot_milan_rat_lags(data_names, axes):

    xpos = [102, 102, 26, 26]  # 75, 75,59, 59
    ymin = 0
    ymax = 1

    for key, n, xpo in zip(axes, data_names, xpos):
        print("rat", n, key, "\n")

        labels = ["", ""]
        data = read_milan_xyy_2col(n)

        if key == "lower left":
            labels = [r"HFA $\rightarrow$ LFP", r"LFP $\rightarrow$ HFA"]

        # if key == "upper right":
        # labels = ["HF amp leads", "LF pha leads"]

        # plot_milan_lags_CMI(axes[key], data, labels)

        if "right" in key:
            plot_milan_lags_MI(axes[key], data, xpo, labels, xlim=200)
            axes[key].axvline(xpo, ymin, ymax, lw=lw + 0.5)
            # axes[key].set_yticks([])
        if "left" in key:
            plot_milan_lags_CMI(axes[key], data, labels)
        if "upper" in key:
            axes[key].set_xticks([])


def read_plot_milan_ross_lags(data, axes):

    xpos = [13, 10, 10, 13]

    for key, d, pos in zip(axes, data, xpos):

        labels = ["", ""]
        print("ross", key)

        axes[key].set_xlim(0, 33)

        if "upper" in key:
            axes[key].set_xticks([])

        if key == "lower left":
            labels = [r"y $\rightarrow$ x", r"x $\rightarrow$ y"]

        # if key == "lower right":
        #    labels = ["y leads", "x leads"]

        if "left" in key:
            plot_milan_lags_CMI(axes[key], d, labels)

        if "right" in key:
            plot_milan_lags_MI(axes[key], d, pos, labels, xlim=33)
            # axes[key].set_yticks([])


def plot_milan_lags_MI(ax, data, xpos, labels, xlim, ymin=0, ymax=1):

    ax.set_xlim(0, xlim)
    ax.axvline(xpos, ymin, ymax, lw=lw + 0.5)
    ax.plot(
        (data[0] - np.min(data[0])) / (np.max(data[0]) - np.min(data[0])),
        color="g",
        ls="-",
        lw=lw,
        label=labels[0],
    )
    ax.plot(
        (data[1] - np.min(data[1])) / (np.max(data[1]) - np.min(data[1])),
        color="r",
        ls="--",
        lw=lw,
        label=labels[1],
    )

    ax.legend(frameon=False, fontsize=fontsize_medium)


def plot_milan_lags_CMI(ax, data, labels):

    ax.plot(
        np.arange(len(data[0]))*5,
        # (data[0]) / (np.max(data[1]) - np.min(data[1])),
        1000 * data[0],
        color="g",
        ls="-",
        lw=lw,
        label=labels[0],
    )
    ax.plot(
        np.arange(len(data[0]))*5,
        # (data[1]) / (np.max(data[1]) - np.min(data[1])),
        1000 * data[1],
        color="r",
        ls="--",
        lw=lw,
        label=labels[1],
    )

    ax.legend(frameon=False, fontsize=fontsize_small - 1)


def prepare_figue_lag_ross(fig, axes):

    fig.text(0.25, 0.9, "CMI", va="center", fontsize=fontsize_big)
    fig.text(0.65, 0.9, "MI", va="center", fontsize=fontsize_big)

    fig.text(0.43, 0.04, "Lag[au]", va="center", fontsize=fontsize_small)
    fig.text(
        0.03, 0.5, "I[nat]", va="center", rotation="vertical", fontsize=fontsize_medium
    )

    # IM-IC LF 8Hz, HF 100 HZ
    # fig.text(0.9, 0.3, 'LF-HF 8-100Hz', va='center', rotation=270, fontsize=10)
    # fig.text(0.9, 0.72, 'LF-HF 8-80Hz', va='center', rotation=270, fontsize=10)

    fig.text(0.905, 0.72, r"$R_1$", va="center", rotation=270, fontsize=fontsize_medium)
    fig.text(0.905, 0.3, r"$R_2$", va="center", rotation=270, fontsize=fontsize_medium)

    fig.text(0.1, 0.9, "x1000", va="center", fontsize=fontsize_small)


def prepare_figue_lag_rat(fig, axes):

    fig.text(0.25, 0.9, "CMI", va="center", fontsize=fontsize_big)
    fig.text(0.65, 0.9, "MI", va="center", fontsize=fontsize_big)

    fig.text(0.43, 0.04, "Lag[su]", va="center", fontsize=fontsize_small)

    # IM-IC LF 8Hz, HF 100 HZ
    fig.text(
        0.9, 0.72, "LF-HF 7-167Hz", va="center", rotation=270, fontsize=fontsize_small
    )
    fig.text(
        0.9, 0.3, "LF-HF 8-100Hz", va="center", rotation=270, fontsize=fontsize_small
    )

    fig.text(0.1, 0.9, "x1000", va="center", fontsize=fontsize_small)
    # fig.text(
    #    0.9, 0.72, "LF-HF 8-80Hz", va="center", rotation=270, fontsize=fontsize_small
    # )


def sort_names(names):
    new_names = ["" for x in range(len(names))]
    for n in names:
        if '5-167_TP' in n:
            new_names[0] = n
        elif '5-167_mi' in n:
            new_names[1] = n
        elif '8-100_TP' in n:
            new_names[2] = n
        elif '8-100_mi' in n:
            new_names[3] = n

    print ('\nNNNNNN', new_names, '\n')
    return new_names


# x10 v1 s1 phAA 12k1LX3ct1 kL4 shi1 su30 v1, 2, 3 variables (Sch-IC, Im-IC, PP-IC)
# x10 v1 s1 AAph 21k1LX3ct1 kL4 shi1 su30


# filename = "../../Desktop/Dropbox/transfer_inormation_prague/plots/fig3_lag/s10v2a8f50TPL3.xyy"
# data = read_milan_xyy_2col(filename)
##plot_milan(data)

# filename = "../../Desktop/Dropbox/transfer_inormation_prague/plots/fig3_lag/r11e31phTPcd4mi2q4.xyy"

# data = read_milan_xyy(filename)
# plot_milan(data)
lw = 2.3
fontsize_small = 14
fontsize_medium = 16
fontsize_big = 18
plt.rc("xtick", labelsize=fontsize_small)
plt.rc("ytick", labelsize=fontsize_small)

"""FIGURE Rat data lag"""
pth1 = "../../Desktop/Dropbox/transfer_inormation_prague/plots/fig3_lag/"
# root = "s10*.xyy"
root = "rat*.xyy"
unsorted_data_names = glob(pth1 + root)
data_names = sort_names(unsorted_data_names)

 #rat_8-100_mi2q8.xyy upper left 
 #rat_5-167_TPL3.xyy upper right 
 #rat_8-100_TPL3.xyy lower left 
 #rat_5-167_mi2q8.xyy lower right 



mossaic_keys = [["upper left", "upper right"], ["lower left", "lower right"]]
fig, axes = plt.subplot_mosaic(
    mossaic_keys,
    # constrained_layout=False,
    # gridspec_kw={"hspace": 0, "wspace": 0},
    # figsize=(15.5, 3.5),
)  #

fig.subplots_adjust(hspace=0)
# fig.subplots_adjust(wspace=0)
read_plot_milan_rat_lags(data_names, axes)
prepare_figue_lag_rat(fig, axes)
plt.savefig(pth1 + "fig3_rat_lags.svg")
plt.savefig(pth1 + "fig3_rat_lags.eps")

"""FIGURE Ross data lag"""
root = "r11*.xyy"
data_names = glob(pth1 + root)
data_resufled = read_milan_xyy_4col(data_names)


mossaic_keys = [["upper left", "upper right"], ["lower left", "lower right"]]
fig, axes = plt.subplot_mosaic(
    mossaic_keys,
    # constrained_layout=False,
    # gridspec_kw={"hspace": 0, "wspace": 0},
)
# fig.subplots_adjust(hspace=0)
plt.rc("xtick", labelsize=fontsize_small)
plt.rc("ytick", labelsize=fontsize_small)
fig.subplots_adjust(hspace=0)
# fig.subplots_adjust(wspace=0)

read_plot_milan_ross_lags(data_resufled, axes)
prepare_figue_lag_ross(fig, axes)
plt.savefig(pth1 + "fig3_ross_lags.svg")
plt.savefig(pth1 + "fig3_ross_lags.eps")
plot_margin = 0.25

x0, x1, y0, y1 = plt.axis()
plt.axis((x0 - plot_margin, x1 + plot_margin, y0 - plot_margin, y1 + plot_margin))
plt.show()
