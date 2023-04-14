#!/usr/bin/env python3

import os
import sys
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter

import numpy as np

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize

from mpl_toolkits.axes_grid1.inset_locator import inset_axes


mpl.rcParams["image.cmap"] = "Spectral_r"
mpl.rcParams["image.origin"] = "lower"
mpl.rcParams["image.interpolation"] = "nearest"
mpl.rcParams["axes.grid"] = True

parser = ArgumentParser(
    description="Plot DD calibration solutions", formatter_class=RawTextHelpFormatter
)
parser.add_argument("sol_file", help="npy solution file")
parser.add_argument("--eff_nr", help="eff solution file", default="eff_nr.npy")
parser.add_argument("--fmin", help="Minimum frequency in MHz", type=float, default=120)
parser.add_argument("--fmax", help="Maximum frequency in MHz", type=float, default=160)
parser.add_argument("--out_dir", help="Output directory", default=".")
parser.add_argument("--out_prefix", help="Prefix to the output filename", default="")
parser.add_argument(
    "--sim",
    help="The solutions are from a simulation so the clusters may not be the default NCP clusters",
    action="store_true",
)

args = parser.parse_args(sys.argv[1:])

stations = [0, 20, 40, 50]

if args.sim:
    clusters = list(range(1, 31))[::2]
    clusters_names = [str(c) for c in clusters]
else:
    clusters = [1, 4, 24, 5, 13, 106, 10, 51, 96, 40, 84, 91, 30, 105]
    clusters_names = [
        "~ NCP",
        "3C61.1",
        "3 deg",
        "4 deg",
        "5 deg",
        "6 deg",
        "7 deg",
        "8 deg",
        "9 deg",
        "10 deg",
        "11 deg",
        "14 deg",
        "Cas A",
        "Cyg A",
    ]

pols = dict(zip(["XX", "YY", "XY", "YX"], [[0, 0], [1, 1], [0, 1], [1, 0]]))
pol_stokes = dict(zip(["I", "V", "U", "Q"], [[0, 0], [1, 1], [0, 1], [1, 0]]))


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


def get_gains(d, cluster, station, eff_nr):
    if cluster > eff_nr.shape[0]:
        return 0
    c_id = slice(int(eff_nr.cumsum()[cluster - 1]), int(eff_nr.cumsum()[cluster]))

    s = d[:, station, :, c_id]
    s = s.transpose(0, 2, 1, 3, 4)
    s = np.concatenate(s)

    return s


def cov2stokes(R):
    I = R[:, :, 0, 0] + R[:, :, 1, 1]
    V = -1j * (R[:, :, 0, 1] - R[:, :, 1, 0])
    Q = R[:, :, 0, 0] - R[:, :, 1, 1]
    U = R[:, :, 0, 1] + R[:, :, 1, 0]

    return I, Q, U, V


def g_mul(g1, g2, c):
    return np.matmul(np.matmul(g1, c), g2.conj().transpose(0, 1, 3, 2))


def get_all_gains(stations, clusters, d, eff_nr):
    gains = dict()
    for s in stations:
        for c in clusters:
            if c > eff_nr.shape[0]:
                print(c, ">", eff_nr.shape[0], "ignoring")
                continue
            gains[(s, c)] = dict()

            g = get_gains(d, c, s, eff_nr)
            gg = g_mul(g, g, np.eye(2))
            I, Q, U, V = cov2stokes(gg)
            gg_stokes = np.stack(
                [np.stack([I, Q], axis=-1), np.stack([U, V], axis=-1)], axis=-1
            )

            for pol_name, (i, j) in pols.items():
                gains[(s, c)][pol_name] = g[:, :, i, j]
                gains[(s, c)]["gg_" + pol_name] = gg[:, :, i, j]
            for pol_name, (i, j) in pol_stokes.items():
                gains[(s, c)]["gg_" + pol_name] = gg_stokes[:, :, i, j]

    return gains


def do_plot(
    freqs,
    gain_data,
    pol_name,
    meta_data,
    out_dir,
    file_name,
    vmin,
    vmax,
    action_fct=np.abs,
    log_norm=True,
):
    fig, axs = plt.subplots(
        ncols=len(stations),
        nrows=len(clusters),
        facecolor="white",
        figsize=(3 * len(stations), 2.5 * len(clusters)),
        sharex=True,
        sharey=False,
        gridspec_kw={"wspace": 0.03, "hspace": 0.03},
    )
    if log_norm:
        norm = LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = Normalize(vmin=vmin, vmax=vmax)

    for i, cluster in enumerate(clusters):
        if (stations[0], cluster) not in gain_data.keys():
            print((stations[0], cluster), "not in", gain_data.keys())
            continue
        for j, station in enumerate(stations):
            g = action_fct(gain_data[(station, cluster)][pol_name])
            g_map = axs[i, j].imshow(
                g,
                aspect="auto",
                extent=[freqs[0], freqs[-1], 0, g.shape[0] - 1],
                norm=norm,
            )
            if j != 0:
                axs[i, j].yaxis.set_ticklabels([])
            if j == 0:
                axs[i, j].set_ylabel("Time slot\n(cluster %s)" % clusters_names[i])
            if i == len(clusters) - 1:
                axs[i, j].set_xlabel("Frequency [MHz]")
            if i == 0:
                axs[i, j].set_title("Station %s" % meta_data["stations"][station])

            if j == len(stations) - 1 and i == 0:
                cbs = ColorbarSetting(ColorbarInnerPosition(height="70%", pad=0.7))
                cbs.add_colorbar(g_map, axs[i, j])

    # fig.tight_layout(pad=0.4)
    fig.savefig(os.path.join(out_dir, "%s.pdf" % file_name), bbox_inches="tight")


def main(args):
    sol_file = args.sol_file
    meta_file_file = sol_file[:-4] + ".npz"
    eff_file = args.eff_nr
    print("Reading %s ..." % sol_file)
    data = np.load(sol_file)

    meta_data = np.load(meta_file_file)
    print(f"{meta_data['stations']} , all {len(meta_data['stations'])} stations")

    eff_nr = np.load(eff_file)

    if eff_nr.sum() != data.shape[4]:
        print("No diffuse emission cluster")
        eff_nr = eff_nr[2:]

    freqs = meta_data["freqs"] * 1e-6
    idx_freqs = (freqs >= args.fmin) & (freqs <= args.fmax)
    freqs = freqs[idx_freqs]
    data = data[:, :, idx_freqs]

    d = data[:, :, :, :, :, 0] + 1j * data[:, :, :, :, :, 1]
    d = d.transpose(0, 1, 2, 4, 3)
    d = d.reshape(*(list(d.shape[:-1]) + [2, 2]))

    print(d.shape, "raw solutions shape")

    print("Computing gains ...")
    gain_data = get_all_gains(stations, clusters, d, eff_nr)

    # print(gain_data.keys(), "computed gains keys")

    print("Starting plotting ...")

    if args.out_prefix:
        args.out_prefix = args.out_prefix + "_"

    for pol_name, (i, j) in pols.items():
        file_name = args.out_prefix + "sol_DD_abs_g_%s" % pol_name
        do_plot(
            freqs,
            gain_data,
            pol_name,
            meta_data,
            args.out_dir,
            file_name,
            1e-1,
            1e1,
            action_fct=np.abs,
            log_norm=True,
        )
        file_name = args.out_prefix + "sol_DD_abs_gg_%s" % pol_name
        do_plot(
            freqs,
            gain_data,
            "gg_" + pol_name,
            meta_data,
            args.out_dir,
            file_name,
            1e-1,
            1e1,
            action_fct=np.abs,
            log_norm=True,
        )
        file_name = args.out_prefix + "sol_DD_phase_g_%s" % pol_name
        do_plot(
            freqs,
            gain_data,
            pol_name,
            meta_data,
            args.out_dir,
            file_name,
            -np.pi,
            np.pi,
            action_fct=np.angle,
            log_norm=False,
        )

    for pol_name, (i, j) in pol_stokes.items():
        file_name = args.out_prefix + "sold_DD_abs_gg_%s" % pol_name
        do_plot(
            freqs,
            gain_data,
            "gg_" + pol_name,
            meta_data,
            args.out_dir,
            file_name,
            5e-1,
            0.5e1,
            action_fct=np.abs,
            log_norm=True,
        )

    print("All done !")


if __name__ == "__main__":
    main(args)
