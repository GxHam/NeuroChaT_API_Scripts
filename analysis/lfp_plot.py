import math
import os
from collections import OrderedDict
import numpy as np

import matplotlib.pyplot as plt

from api_utils import make_dir_if_not_exists
from neurochat.nc_lfp import NLfp
from neurochat.nc_utils import butter_filter


def plot_long_lfp(
        lfp, out_name, nsamples=None, offset=0,
        nsplits=3, figsize=(32, 4), ylim=(-0.4, 0.4)):
    """
    Create a figure to display a long LFP signal in nsplits rows.

    Args:
        lfp (NLfp): The lfp signal to plot.
        out_name (str): The name of the output, including directory.
        nsamples (int, optional): Defaults to all samples.
            The number of samples to plot. 
        offset (int, optional): Defaults to 0.
            The number of samples into the lfp to start plotting at.
        nsplits (int, optional): The number of rows in the resulting figure.
        figsize (tuple of int, optional): Defaults to (32, 4)
        ylim (tuple of float, optional): Defaults to (-0.4, 0.4)

    Returns:
        None

    """
    if nsamples is None:
        nsamples = lfp.get_total_samples()

    fig, axes = plt.subplots(nsplits, 1, figsize=figsize)
    for i in range(nsplits):
        start = int(offset + i * (nsamples // nsplits))
        end = int(offset + (i + 1) * (nsamples // nsplits))
        if nsplits == 1:
            ax = axes
        else:
            ax = axes[i]
        ax.plot(
            lfp.get_timestamp()[start:end],
            lfp.get_samples()[start:end], color='k')
        ax.set_ylim(ylim)
    plt.tight_layout()
    fig.savefig(out_name, dpi=400)
    plt.close(fig)


def plot_lfp(
        out_dir, lfp_odict, segment_length=150, in_range=None, dpi=50):
    """
    Create a number of figures to display lfp signal on multiple channels.

    There will be one figure for each split of thetotal length 
    into segment_length, and a row for each value in lfp_odict.

    It is assumed that the input lfps are prefiltered if filtering is required.

    Args:
        out_dir (str): The name of the file to plot to, including dir.
        lfp_odict (OrderedDict): Keys as channels and NLfp objects as vals.
        segment_length (float): Time in seconds of LFP to plot in each figure.
        in_range (tuple(int, int), optional): Time(s) of LFP to plot overall.
            Defaults to None.
        dpi (int, optional): Resulting plot dpi.

    Returns:
        None

    """
    if in_range is None:
        in_range = (0, max([lfp.get_duration() for lfp in lfp_odict.values()]))
    y_axis_max = max([max(lfp.get_samples()) for lfp in lfp_odict.values()])
    y_axis_min = min([min(lfp.get_samples()) for lfp in lfp_odict.values()])
    for split in np.arange(in_range[0], in_range[1], segment_length):
        fig, axes = plt.subplots(
            nrows=len(lfp_odict),
            figsize=(20, len(lfp_odict) * 2))
        a = np.round(split, 2)
        b = np.round(min(split + segment_length, in_range[1]), 2)
        out_name = os.path.join(
            out_dir, "{}s_to_{}s.png".format(a, b))
        for i, (key, lfp) in enumerate(lfp_odict.items()):
            convert = lfp.get_sampling_rate()
            c_start, c_end = math.floor(a * convert), math.floor(b * convert)
            lfp_sample = lfp.get_samples()[c_start:c_end]
            x_pos = lfp.get_timestamp()[c_start:c_end]
            axes[i].plot(x_pos, lfp_sample, c="k")
            axes[i].text(
                0.03, 1, "Channel " + key,
                transform=axes[i].transAxes, c="k")
            axes[i].set_ylim(y_axis_min, y_axis_max)
            axes[i].set_xlim(a, b)
        print("Saving result to {}".format(out_name))
        fig.savefig(out_name, dpi=dpi)
        plt.close("all")


def plot_sample_of_signal(
        load_loc, out_dir=None, name=None, offseta=0, length=50,
        filt_params=(False, None, None)):
    """
    Plot a small filtered sample of the LFP signal in the given band.

    offseta and length are times
    """
    in_dir = os.path.dirname(load_loc)
    lfp = NLfp()
    lfp.load(load_loc)

    if out_dir is None:
        out_loc = "nc_signal"
        out_dir = os.path.join(in_dir, out_loc)

    if name is None:
        name = "full_signal_filt.png"

    make_dir_if_not_exists(out_dir)
    out_name = os.path.join(out_dir, name)
    fs = lfp.get_sampling_rate()
    filt, lower, upper = filt_params
    lfp_samples = lfp.get_samples()
    if filt:
        lfp_samples = butter_filter(
            lfp_samples, fs, 10, lower, upper, 'bandpass')
    plot_long_lfp(
        lfp_samples, out_name, nsplits=1, ylim=(-0.3, 0.3), figsize=(20, 8),
        offset=lfp.get_sampling_rate() * offseta,
        nsamples=lfp.get_sampling_rate() * length)


def plot_coherence(f, Cxy, name=None, dpi=100):
    fig, ax = plt.subplots()
    # ax.semilogy(f, Cxy)
    ax.plot(f, Cxy, c="k")
    ax.set_xlabel('frequency [Hz]')
    ax.set_ylabel('Coherence')
    ax.set_xticks(np.arange(0, f.max(), 10))
    ax.set_ylim(0, 1)
    if name is None:
        plt.show()
    else:
        fig.savefig(name, dpi=dpi)
