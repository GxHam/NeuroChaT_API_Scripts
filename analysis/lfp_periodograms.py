import os

from api_utils import read_cfg, parse_args, setup_logging, make_dir_if_not_exists, make_path_if_not_exists
from lfp_odict import LfpODict
import neurochat.nc_plot as nc_plot
from neurochat.nc_utils import get_all_files_in_dir

import matplotlib.pyplot as plt
from lfp_plot import plot_lfp

import api_plot_org as plot_org


def main(fname, analysis_flags, o_main_dir=None):
    '''
    Parameters
    ----------
    fname : str
        filenames to be analysed

    analysis_flags : bool, optional. Defaults to True.
        Sets analysis to be used.
        0 - plot periodograms and ptrs in seperate plots for each tetrode
        1 - plot graphs from all tetrodes in 1 .png

    o_main_dir: dir, optional. Defaults to None.
        None - Saves plots in a LFP folder where .eeg was found
        Else - Saves plots in a LFP folder of given drive
    '''

    # Setup drive and LFPs to extract
    chans = [i for i in range(1, 17)]
    lfp_odict = LfpODict(
        fname, channels=chans, filt_params=(True, 1.5, 90))
    if o_main_dir is None:
        o_dir = os.path.join(
            os.path.dirname(fname), "!LFP")
    else:
        o_dir = os.path.join(o_main_dir, "!LFP")
    make_dir_if_not_exists(o_dir)

    if analysis_flags[0]:   # Plot periodograms and ptr for each tetrode seperately
        # Plot periodogram for each eeg
        for i, (key, lfp) in enumerate(lfp_odict.get_signal().items()):
            graph_data = lfp.spectrum(
                ptype='psd', prefilt=False,
                db=False, tr=False,
                filtset=[10, 1.0, 40, 'bandpass'])
            fig = nc_plot.lfp_spectrum(graph_data)
            plt.ylim(0, 0.01)
            plt.xlim(0, 40)
            out_name = os.path.join(o_dir, "p", key + "p.png")
            make_path_if_not_exists(out_name)
            fig.savefig(out_name)
            plt.close()

            graph_data = lfp.spectrum(
                ptype='psd', prefilt=False,
                db=True, tr=True,
                filtset=[10, 1.0, 40, 'bandpass'])
            fig = nc_plot.lfp_spectrum_tr(graph_data)
            # plt.ylim(0, 0.01)
            # plt.xlim(0, 40)
            out_name = os.path.join(o_dir, "ptr", key + "ptr.png")
            make_path_if_not_exists(out_name)
            print("Saving result to {}".format(out_name))
            fig.savefig(out_name)
            plt.close()

        plot_lfp(o_dir, lfp_odict.get_filt_signal(), segment_length=60)

    if analysis_flags[1]:   # Complie graphs per session in a single .png
        # Region info for eeg
        names = ["CLA"] * 8 + ["ACC"] * 4 + ["RSC"] * 4

        # Setup summary grid
        rows, cols = [4, 4]
        gf = plot_org.GridFig(rows, cols, wspace=0.3,
                              hspace=0.3, tight_layout=False)

        # Plot summary periodogram
        for i, (key, lfp) in enumerate(lfp_odict.get_signal().items()):
            graph_data = lfp.spectrum(
                ptype='psd', prefilt=False,
                db=False, tr=False,
                filtset=[10, 1.0, 40, 'bandpass'])
            ax = gf.get_next(along_rows=False)
            nc_plot.lfp_spectrum(graph_data, ax)
            plt.ylim(0, 0.015)
            # plt.xlim(0, 40)
            if i % 4 == 0:
                ax.text(0.49, 1.08, names[i], fontsize=20,
                        horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        gf.fig.suptitle(
            (fname.split("\\")[-1][4:] + " Periodogram"), fontsize=30)
        out_name = os.path.join(
            o_dir, "Sum_p", fname.split("\\")[-1] + "_p_sum.png")
        make_path_if_not_exists(out_name)
        print("Saving result to {}".format(out_name))
        gf.fig.savefig(out_name)
        plt.close()

        # Plot summary periodogram tr
        gf = plot_org.GridFig(rows, cols, wspace=0.5, hspace=0.5)
        for i, (key, lfp) in enumerate(lfp_odict.get_signal().items()):
            graph_data = lfp.spectrum(
                ptype='psd', prefilt=True,
                db=True, tr=True,
                filtset=[10, 1.0, 40, 'bandpass'])
            ax = gf.get_next(along_rows=False)
            nc_plot.lfp_spectrum_tr(graph_data, ax)
            plt.ylim(0, 40)
            # plt.xlim(0, 40)
            if i % 4 == 0:
                ax.text(0.49, 1.08, names[i], fontsize=20,
                        horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

        gf.fig.suptitle(
            (fname.split("\\")[-1][4:] + " Time Resolved Periodogram"), fontsize=30)
        out_name = os.path.join(
            o_dir, "Sum_ptr", fname.split("\\")[-1] + "_ptr_sum.png")
        make_path_if_not_exists(out_name)
        print("Saving result to {}".format(out_name))
        gf.fig.savefig(out_name)
        plt.close()


if __name__ == "__main__":
    in_dir = r"F:\Ham Data\A10_CAR-SA2"

    filenames = get_all_files_in_dir(
        in_dir, ext=".eeg", recursive=True,
        verbose=True, re_filter="Pre")

    filenames = [fname[:-4] for fname in filenames]
    if len(filenames) == 0:
        print("No set files found for analysis!")
        exit(-1)

    # Description of analysis_flags
    # 0 - plot periodograms and ptrs in seperate plots for each tetrode
    # 1 - plot graphs from all tetrodes in 1 .png
    # 2 -
    # 3 -

    analysis_flags = [0, 1, 0, 0]
    for fname in filenames:
        main(fname, analysis_flags, in_dir)
