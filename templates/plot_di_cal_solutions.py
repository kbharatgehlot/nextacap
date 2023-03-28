#!/usr/bin/env python3

import os
import sys
import collections
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter

# from statsmodels import robust
from astropy.stats import median_absolute_deviation as mad

import numpy as np

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from mpl_toolkits.axes_grid1.inset_locator import inset_axes


mpl.rcParams["image.cmap"] = "Spectral_r"
mpl.rcParams["image.origin"] = "lower"
mpl.rcParams["image.interpolation"] = "nearest"
mpl.rcParams["axes.grid"] = True

parser = ArgumentParser(
    description="Plot DI calibration solutions", formatter_class=RawTextHelpFormatter
)
parser.add_argument("sol_file", help="npy solution file")
parser.add_argument("--fmin", help="Minimum frequency in MHz", type=float, default=120)
parser.add_argument("--fmax", help="Maximum frequency in MHz", type=float, default=160)
parser.add_argument("--cluster", help="Cluster index", default=0, type=int)
parser.add_argument("--out_dir", help="Output directory", default=".")
parser.add_argument("--out_prefix", help="Prefix to the output filename", default="")

station_pairs = [(0, 1), (5, 8), (15, 18), (36, 39), (52, 55), (55, 59)]


class ColorbarInnerPosition(object):
    def __init__(
        self,
        orientation="horizontal",
        width="5%",
        height="50%",
        location=1,
        pad=0.5,
        tick_position=None,
    ):
        """
        width, height: inch if number, percentage of parent axes if string (like '5%')
        pad: points
        location are :
        'upper right' : 1,
        'upper left' : 2,
        'lower left' : 3,
        'lower right' : 4,
        'right' : 5,
        'center left' : 6,
        'center right' : 7,
        'lower center' : 8,
        'upper center' : 9,
        'center' : 10,
        """
        self.orientation = orientation
        if orientation == "vertical":
            self.width = width
            self.height = height
            if tick_position is None:
                tick_position = "left"
        else:
            self.width = height
            self.height = width
            if tick_position is None:
                tick_position = "bottom"
        self.location = location
        self.pad = pad
        self.tick_position = tick_position

    def get_cb_axes(self, ax):
        cax = inset_axes(
            ax,
            width=self.width,
            height=self.height,
            loc=self.location,
            borderpad=self.pad,
        )
        return cax

    def post_creation(self, colorbar):
        if self.orientation == "vertical":
            if self.tick_position == "left":
                colorbar.ax.yaxis.set_ticks_position(self.tick_position)
                colorbar.ax.yaxis.set_label_position(self.tick_position)
        else:
            if self.tick_position == "top":
                colorbar.ax.xaxis.set_ticks_position(self.tick_position)
                colorbar.ax.xaxis.set_label_position(self.tick_position)

    def get_orientation(self):
        return self.orientation


class ColorbarSetting(object):
    def __init__(
        self, cb_position, ticks_locator=None, ticks_formatter=None, cmap="jet"
    ):
        self.cb_position = cb_position
        self.ticks_locator = ticks_locator
        self.ticks_formatter = ticks_formatter
        self.cmap = cmap

    def add_colorbar(self, mappable, ax):
        fig = ax.get_figure()
        cb = fig.colorbar(
            mappable,
            ticks=self.ticks_locator,
            format=self.ticks_formatter,
            orientation=self.cb_position.get_orientation(),
            cax=self.cb_position.get_cb_axes(ax),
        )
        self.cb_position.post_creation(cb)
        if not hasattr(fig, "_plotutils_colorbars"):
            fig._plotutils_colorbars = dict()
        fig._plotutils_colorbars[ax] = cb
        return cb

    def get_cmap(self):
        return self.cmap


def g_mul(g1, g2, c):
    return np.matmul(np.matmul(g1, c), g2.conj().transpose(0, 1, 3, 2))


def cov2stokes(R):
    I = R[:, :, :, 0, 0] + R[:, :, :, 1, 1]
    V = -1j * (R[:, :, :, 0, 1] - R[:, :, :, 1, 0])
    Q = R[:, :, :, 0, 0] - R[:, :, :, 1, 1]
    U = R[:, :, :, 0, 1] + R[:, :, :, 1, 0]

    return I, Q, U, V


def plot_abs_and_phase(freqs, d, pols, out_dir, prefix):
    fig, axs = plt.subplots(
        ncols=4,
        nrows=len(station_pairs),
        figsize=(12, 2.5 * len(station_pairs)),
        sharex=True,
        sharey=False,
    )
    fig2, axs2 = plt.subplots(
        ncols=4,
        nrows=len(station_pairs),
        figsize=(12, 2.5 * len(station_pairs)),
        sharex=True,
        sharey=False,
    )

    c = np.array([[1, 0], [0, 1]])
    myestimator = np.abs(g_mul(d[:, 0], d[:, 1], c)[:, :, 0, 0])
    print(np.average(np.log10(myestimator)), np.std(np.log10(myestimator)))

    myvmin = 10 ** (
        np.int(
            np.average(np.log10(myestimator)) - 3 * np.std(np.log10(myestimator)) - 0.5
        )
    )
    myvmax = 10 ** (
        np.int(
            np.average(np.log10(myestimator)) + 3 * np.std(np.log10(myestimator)) + 0.5
        )
    )
    if myvmax == myvmin:
        myvmax = 10 * myvmin
    for i, (s1, s2) in enumerate(station_pairs):
        print(d.shape)
        g1 = d[:, s1]
        g2 = d[:, s2]

        for j, (name, (k, l)) in enumerate(pols.items()):
            extent = [freqs[0], freqs[-1], 0, g1.shape[0] - 1]
            gg = g_mul(g1, g2, c)[:, :, k, l]

            g_map = axs[i, j].imshow(
                abs(gg),
                aspect="auto",
                norm=LogNorm(),
                vmin=myvmin,
                vmax=myvmax,
                extent=extent,
            )
            p_map = axs2[i, j].imshow(
                np.angle(gg), aspect="auto", vmin=-np.pi, vmax=np.pi, extent=extent
            )

            for ax in [axs, axs2]:
                if j != 0:
                    ax[i, j].yaxis.set_ticklabels([])
                if j == 0:
                    ax[i, j].set_ylabel("Time slot\n(stations %s - %s)" % (s1, s2))
                if i == len(station_pairs) - 1:
                    ax[i, j].set_xlabel("Frequency [MHz]")
                if i == 0:
                    ax[i, j].set_title("g1g2*_%s" % name)
                if j == len(pols) - 1 and i == 0:
                    cbs = ColorbarSetting(ColorbarInnerPosition(height="70%", pad=0.7))
                    cbs.add_colorbar(g_map, axs[i, j])
                    cbs = ColorbarSetting(ColorbarInnerPosition(height="70%", pad=0.7))
                    cbs.add_colorbar(p_map, axs2[i, j])

    fig.suptitle("Sagecal solutions abs(g1g2*)", va="bottom", y=1)
    fig2.suptitle("Sagecal solutions phase(g1g2*)", va="bottom", y=1)

    fig.tight_layout(pad=0.4)
    fig2.tight_layout(pad=0.4)

    fig.savefig(os.path.join(out_dir, "%s_abs_gg.pdf" % prefix))
    fig2.savefig(os.path.join(out_dir, "%s_phase_gg.pdf" % prefix))


def plot_sol_mean_and_std(freqs, gg, pols, out_dir, avg_time=True, prefix="sol_DI"):
    if avg_time:
        avg_col = "time"
        avg_col_idx = 1
        other_col = "Frequency [MHz]"
        extent = [freqs[0], freqs[-1], 0, gg.shape[0] - 1]
    else:
        avg_col = "frequency"
        avg_col_idx = 2
        other_col = "Time index"
        extent = [0, gg.shape[1] - 1, 0, gg.shape[0] - 1]

    fig, axs = plt.subplots(
        ncols=4, nrows=3, figsize=(12, 7.5), sharex=True, sharey=True
    )
    myestimator = np.median(np.abs(gg[:, :, :, 0, 0]), axis=avg_col_idx)
    myvmin = 10 ** (
        np.int(
            np.average(np.log10(myestimator)) - 3 * np.std(np.log10(myestimator)) - 0.5
        )
    )
    myvmax = 10 ** (
        np.int(
            np.average(np.log10(myestimator)) + 3 * np.std(np.log10(myestimator)) + 0.5
        )
    )

    for j, (name, (k, l)) in enumerate(pols.items()):
        mean_map = axs[0, j].imshow(
            np.median(np.abs(gg[:, :, :, k, l]), axis=avg_col_idx),
            aspect="auto",
            norm=LogNorm(),
            vmin=myvmin,
            vmax=myvmax,
            extent=extent,
        )
        rms_map = axs[1, j].imshow(
            mad(np.real(gg[:, :, :, k, l]), axis=avg_col_idx),
            aspect="auto",
            norm=LogNorm(),
            vmin=myvmin,
            vmax=myvmax,
            extent=extent,
        )

        gg_diff = np.diff(gg[:, :, :, k, l], axis=2)
        rms_diff_map = axs[2, j].imshow(
            0.5 * mad(np.real(gg_diff), axis=avg_col_idx),
            aspect="auto",
            norm=LogNorm(),
            vmin=myvmin,
            vmax=myvmax,
            extent=extent,
        )

        if j != 0:
            axs[0, j].yaxis.set_ticklabels([])
            axs[1, j].yaxis.set_ticklabels([])

        axs[0, 0].set_ylabel("Mean along %s\n\nStations" % avg_col)
        axs[1, 0].set_ylabel("Std along %s\n\nStations" % avg_col)
        axs[2, 0].set_ylabel("Std along %s\nof the SB diff gains\n\nStations" % avg_col)
        axs[2, j].set_xlabel(other_col)
        axs[0, j].set_title("gg*_%s" % name)

        if j == len(pols) - 1:
            cbs = ColorbarSetting(ColorbarInnerPosition(height="70%", pad=0.7))
            cbs.add_colorbar(mean_map, axs[0, j])
            cbs = ColorbarSetting(ColorbarInnerPosition(height="70%", pad=0.7))
            cbs.add_colorbar(rms_map, axs[1, j])
            cbs = ColorbarSetting(ColorbarInnerPosition(height="70%", pad=0.7))
            cbs.add_colorbar(rms_diff_map, axs[2, j])

    fig.tight_layout(pad=0.4)
    fig.savefig(os.path.join(out_dir, "%s_mean_std_%s_gg.pdf" % (prefix, avg_col)))


def main():
    args = parser.parse_args(sys.argv[1:])
    sol_file = args.sol_file
    meta_file_file = sol_file[:-4] + ".npz"

    print("Reading %s ..." % sol_file)
    data = np.load(sol_file)
    meta_data = np.load(meta_file_file)

    # print("# stations", len(meta_data["stations"]))

    freqs = meta_data["freqs"] * 1e-6
    idx_freqs = (freqs >= args.fmin) & (freqs <= args.fmax)
    freqs = freqs[idx_freqs]
    data = data[:, :, idx_freqs]

    d = data[:, :, :, :, args.cluster, 0] + 1j * data[:, :, :, :, args.cluster, 1]
    d = d.reshape(*(list(d.shape[:-1]) + [2, 2]))

    print("Computing gains ...")
    c = np.array([[1, 0], [0, 1]])
    gg = np.array([g_mul(d[:, i], d[:, i], c) for i in range(d.shape[1])])
    pols = collections.OrderedDict(
        zip(["XX", "YY", "XY", "YX"], [[0, 0], [1, 1], [0, 1], [1, 0]])
    )

    I, Q, U, V = cov2stokes(gg)
    gg_stokes = np.stack(
        [np.stack([I, Q], axis=-1), np.stack([U, V], axis=-1)], axis=-1
    )
    pol_stokes = collections.OrderedDict(
        zip(["I", "V", "U", "Q"], [[0, 0], [1, 1], [0, 1], [1, 0]])
    )

    print("Starting plotting ...")

    if args.out_prefix:
        args.out_prefix = args.out_prefix + "_"

    plot_abs_and_phase(freqs, d, pols, args.out_dir, prefix=args.out_prefix + "sol_DI")

    plot_sol_mean_and_std(
        freqs, gg, pols, args.out_dir, avg_time=True, prefix=args.out_prefix + "sol_DI"
    )
    plot_sol_mean_and_std(
        freqs, gg, pols, args.out_dir, avg_time=False, prefix=args.out_prefix + "sol_DI"
    )

    plot_sol_mean_and_std(
        freqs,
        gg_stokes,
        pol_stokes,
        args.out_dir,
        avg_time=True,
        prefix=args.out_prefix + "sol_stokes_DI",
    )
    plot_sol_mean_and_std(
        freqs,
        gg_stokes,
        pol_stokes,
        args.out_dir,
        avg_time=False,
        prefix=args.out_prefix + "sol_stokes_DI",
    )

    print("All done!")


if __name__ == "__main__":
    main()
