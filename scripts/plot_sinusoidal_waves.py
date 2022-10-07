import pandas as pd
import mne

from mne_connectivity import SpectralConnectivity
from itertools import product
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

sns.set_style("white")
# matplotlib inline


def sin_waves(num_waves, amp=1):

    x = np.arange(0, num_waves * 2 * np.pi)
    y = amp * np.sin(x)

    return x, y


# x, y = sin_waves()
# fig, ax = plt.subplots()
# ax.plot(x, y)


# We can generate these sine wave parameters, then stitch them together
amplitude_values = [1, 3, 10, 15]
phase_values = [0, 0.25, 0.33, 0.5]
freq = 2
signal_vals = list(product(amplitude_values, phase_values))
amps, phases = zip(*signal_vals)

# We'll also define some noise levels to see how this affects results
noise_levels = [0, 2, 4, 8]

# Now define how long these signals will be
t_stop = 50
time = np.arange(0, t_stop, 0.01)

# We're storing everything in dataframes, so create some indices
ix_amp = pd.Index(amps, name="amp")
ix_phase = pd.Index(phases, name="phase")


def create_sin_signals():
    # Create all our signals
    signals = []
    for noise_level in noise_levels:
        sig_ = np.array(
            [
                amp * np.sin(freq * 2 * np.pi * time + 2 * np.pi * phase)
                for amp, phase in signal_vals
            ]
        )
        noise = noise_level * np.random.randn(*sig_.shape)
        sig_ += noise
        ix_noise = pd.Index([noise_level] * sig_.shape[0], name="noise_level")
        ix_multi = pd.MultiIndex.from_arrays([ix_amp, ix_phase, ix_noise])
        signals.append(pd.DataFrame(sig_, index=ix_multi))
    signals = pd.concat(signals, 0)
    signals.columns.name = "time"
    return signals, noise


def compute_connectivity(signals):
    con_all = []
    for ix_noise, sig in signals.groupby(level="noise_level"):
        # Setting up output indices
        this_noise_level = sig.index.get_level_values("noise_level").unique()[0]
        ix_ref = np.where(sig.eval("amp==3 and phase==0"))[0][0]
        ix_time = pd.Index(range(sig.shape[0]), name="time")
        ix_cc = pd.Index(["cc"] * sig.shape[0], name="con")
        ix_coh = pd.Index(["coh"] * sig.shape[0], name="con")

        # Calculating correlation is easy with pandas
        cc = sig.T.corr().astype(float).iloc[:, ix_ref]
        cc.name = None
        cc = pd.DataFrame(cc)

        # We'll use MNE for coherenece
        indices = (np.arange(sig.shape[0]), ix_ref.repeat(sig.shape[0]))
        # con, freqs, times, epochs, tapers = SpectralConnectivity(
        #    sig.values[None, :, :], sfreq=freq, indices=indices)
        # con_mn = con.mean(-1)
        # con_mn = pd.DataFrame(con_mn, index=cc.index)

        # Final prep
        # con_mn = con_mn.set_index(ix_coh, append=True)
        cc = cc.set_index(ix_cc, append=True)
        # con_all += ([con_mn, cc])
    # con_all = pd.concat(con_all, axis=0).squeeze().unstack('noise_level')

    return ix_ref, ix_cc, cc  # con_all


def plot_coherence_corr_sig(con_all, noise):

    f, axs = plt.subplots(2, 2, figsize=(15, 10))
    for ax, (noise, vals) in zip(axs.ravel(), con_all.iteritems()):
        ax = vals.unstack("con").plot(ax=ax)
        ax.set_title("Noise level: {0}".format(noise))
        ax.set_ylim([-1.1, 1.1])


def plot_coherence_vs_corr():

    plt_df = con_all.stack("noise_level").unstack("con").reset_index("noise_level")
    ax = plt_df.plot(
        "cc",
        "coh",
        c="noise_level",
        kind="scatter",
        cmap=plt.cm.Reds,
        figsize=(10, 5),
        alpha=0.5,
        s=50,
    )
    ax.set_title("CC vs Coherence")


def plot_signals(n_combs, sig_combinations, noise_level, cc):
    fig, axs = plt.subplots(int(n_combs / 4), 4, figsize=(15, n_combs / 3 * 5))
    plt_sig = signals.xs(noise_level, level="noise_level")
    n_combs = len(sig_combinations)

    for (comp_a, comp_b), ax in zip(sig_combinations, axs.ravel()):
        plt_sig.iloc[[comp_a, comp_b]].T.head(250).plot(ax=ax, legend=None)
        # ax.set_title("CC: {0}".format(c))


signals, noise = create_sin_signals()
ix_ref, ix_cc, cc = compute_connectivity(signals)
# plot_coherence_corr_sig(con_all, noise)

# Set up a dataframe for plotting
noise_level = noise_levels[1]

# Define 16 signals to plot
sig_combinations = list(product([ix_ref], range(16)))
n_combs = len(sig_combinations)
plot_signals(n_combs, sig_combinations, noise_level, cc)


plt.show()
