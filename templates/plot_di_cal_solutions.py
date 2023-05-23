#!/usr/bin/env python

from __future__ import print_function
from __future__ import absolute_import

import os
import re
import sys
import collections
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter

import astropy.stats as astats

import numpy as np

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from mpl_toolkits.axes_grid1.inset_locator import inset_axes


mpl.rcParams['image.cmap'] = 'Spectral_r'
mpl.rcParams['image.origin'] = 'lower'
mpl.rcParams['image.interpolation'] = 'nearest'
mpl.rcParams['axes.grid'] = True

parser = ArgumentParser(description="Plot DI calibration solutions statistics",
                        formatter_class=RawTextHelpFormatter)
parser.add_argument('sol_file', help='npy solution file')
parser.add_argument('--fmin', help='Minimum frequency in MHz', type=float, default=120)
parser.add_argument('--fmax', help='Maximum frequency in MHz', type=float, default=160)
parser.add_argument('--cluster', help='Cluster index', default=0, type=int)
parser.add_argument('--out_dir', help='Output directory', default='.')
parser.add_argument('--out_prefix', help='Prefix to the output filename', default='')

station_pairs = [(0, 1), (5, 8), (15, 18), (36, 39), (52, 55), (55, 59)]


class ColorbarInnerPosition(object):

    def __init__(self, orientation="horizontal", width="5%", height="50%", location=1, pad=0.5,
                 tick_position=None):
        '''
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
        '''
        self.orientation = orientation
        if orientation == 'vertical':
            self.width = width
            self.height = height
            if tick_position is None:
                tick_position = 'left'
        else:
            self.width = height
            self.height = width
            if tick_position is None:
                tick_position = 'bottom'
        self.location = location
        self.pad = pad
        self.tick_position = tick_position

    def get_cb_axes(self, ax):
        cax = inset_axes(ax, width=self.width, height=self.height, loc=self.location, borderpad=self.pad)
        return cax

    def post_creation(self, colorbar):
        if self.orientation == 'vertical':
            if self.tick_position == 'left':
                colorbar.ax.yaxis.set_ticks_position(self.tick_position)
                colorbar.ax.yaxis.set_label_position(self.tick_position)
        else:
            if self.tick_position == 'top':
                colorbar.ax.xaxis.set_ticks_position(self.tick_position)
                colorbar.ax.xaxis.set_label_position(self.tick_position)

    def get_orientation(self):
        return self.orientation


class ColorbarSetting(object):

    def __init__(self, cb_position, ticks_locator=None, ticks_formatter=None, cmap='jet'):
        self.cb_position = cb_position
        self.ticks_locator = ticks_locator
        self.ticks_formatter = ticks_formatter
        self.cmap = cmap

    def add_colorbar(self, mappable, ax):
        fig = ax.get_figure()
        cb = fig.colorbar(mappable, ticks=self.ticks_locator, format=self.ticks_formatter,
                          orientation=self.cb_position.get_orientation(), cax=self.cb_position.get_cb_axes(ax))
        self.cb_position.post_creation(cb)
        if not hasattr(fig, '_plotutils_colorbars'):
            fig._plotutils_colorbars = dict()
        fig._plotutils_colorbars[ax] = cb
        return cb

    def get_cmap(self):
        return self.cmap


def g_mul(g1, g2, c):
    return np.matmul(np.matmul(g1, c), g2.conj().transpose(0, 1, 3, 2))


def cov2stokes(R):
    I = 0.5 * (R[:, :, :, 0, 0] + R[:, :, :, 1, 1])
    V = 0.5 * (-1j * (R[:, :, :, 0, 1] - R[:, :, :, 1, 0]))
    Q = 0.5 * (R[:, :, :, 0, 0] - R[:, :, :, 1, 1])
    U = 0.5 * (R[:, :, :, 0, 1] + R[:, :, :, 1, 0])

    return I, Q, U, V


def avg_gg(a, n=2):
    end = n * int(a.shape[1] / n)
    return np.mean(a[:, :end, :].reshape((a.shape[0], -1, n, a.shape[2])), axis=2)


def plot_abs_and_phase(freqs, d, pols, out_dir, prefix):
    fig, axs = plt.subplots(ncols=4, nrows=len(station_pairs), figsize=(12, 2.5 * len(station_pairs)),
                            sharex=True, sharey=False)
    fig2, axs2 = plt.subplots(ncols=4, nrows=len(station_pairs), figsize=(12, 2.5 * len(station_pairs)),
                              sharex=True, sharey=False)

    for i, (s1, s2) in enumerate(station_pairs):
        g1 = d[:, s1]
        g2 = d[:, s2]
        c = np.array([[1, 0], [0, 1]])

        for j, (name, (k, l)) in enumerate(pols.items()):
            extent = [freqs[0], freqs[-1], 0, g1.shape[0] - 1]
            gg = g_mul(g1, g2, c)[:, :, k, l]

            if i == 0 and j == 0:
                vmax = np.median(abs(gg)) + 10 * astats.mad_std(abs(gg).flatten())
                vmin = 1e-1 * vmax

            g_map = axs[i, j].imshow(abs(gg), aspect='auto', norm=LogNorm(vmin=vmin, vmax=vmax), extent=extent)
            p_map = axs2[i, j].imshow(np.angle(gg), aspect='auto', vmin=-np.pi, vmax=np.pi, extent=extent)

            axs[i, j].text(0.05, 0.95, 'Median: %.4g' % np.median(abs(gg)), transform=axs[i,j].transAxes, va='top', ha='left')
            axs2[i, j].text(0.05, 0.95, 'Median: %.4g' % np.median(np.angle(gg)), transform=axs2[i,j].transAxes, va='top', ha='left')
            
            for ax in [axs, axs2]:
                if j != 0:
                    ax[i, j].yaxis.set_ticklabels([])
                if j == 0:
                    ax[i, j].set_ylabel('Time slot\n(stations %s - %s)' % (s1, s2))
                if i == len(station_pairs) - 1:
                    ax[i, j].set_xlabel('Frequency [MHz]')
                if i == 0:
                    ax[i, j].set_title('g1g2*_%s' % name)
                if j == len(pols) - 1 and i == 0:
                    cbs = ColorbarSetting(ColorbarInnerPosition(height='70%', pad=0.7))
                    cbs.add_colorbar(g_map, axs[i, j])
                    cbs = ColorbarSetting(ColorbarInnerPosition(height='70%', pad=0.7))
                    cbs.add_colorbar(p_map, axs2[i, j])

    fig.suptitle('Sagecal solutions abs(g1g2*)', va='bottom', y=1)
    fig2.suptitle('Sagecal solutions phase(g1g2*)', va='bottom', y=1)

    fig.tight_layout(pad=0.4)
    fig2.tight_layout(pad=0.4)

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    fig.savefig(os.path.join(out_dir, '%s_abs_gg.pdf' % prefix))
    fig2.savefig(os.path.join(out_dir, '%s_phase_gg.pdf' % prefix))


def plot_sol_mean_and_std(freqs, gg, pols, out_dir, avg_type='stations', n_time_avg=40,
                          max_station=62, prefix='sol_DI'):
    if avg_type == 'stations':
        avg_col_idx = 0
        other_col = 'Frequency [MHz]'
        extent = [freqs[0], freqs[-1], 0, gg.shape[1] - 1]
    elif avg_type == 'frequency':
        avg_col_idx = 2
        other_col = 'Stations'
        extent = [0, gg.shape[0] - 1, 0, gg.shape[1] - 1]

    fig, axs = plt.subplots(ncols=4, nrows=4, figsize=(12, 10), sharex=True, sharey=False)

    for j, (name, (k, l)) in enumerate(pols.items()):
        g_m = np.median(np.abs(gg[:max_station, :, :, k, l]), axis=avg_col_idx)
        g_d_t = astats.mad_std(np.real(np.diff(gg[:max_station, :, :, k, l], axis=1)), axis=avg_col_idx) * np.sqrt(0.5)
        g_d_f = astats.mad_std(np.real(np.diff(gg[:max_station, :, :, k, l], axis=2)), axis=avg_col_idx) * np.sqrt(0.5)
        g_d_f_avg = astats.mad_std(np.real(avg_gg(np.diff(gg[:max_station, :, :, k, l], axis=2), n_time_avg)),
                               axis=avg_col_idx) * np.sqrt(0.5 * n_time_avg)

        if j == 0:
            vmax = np.median(g_m) + 10 * astats.mad_std(g_m.flatten())
            vmin = 1e-3 * vmax

        for ax, a in zip(axs[:, j], [g_m, g_d_t, g_d_f, g_d_f_avg]):
            if avg_type == 'frequency':
                a = a.T

            im_map = ax.imshow(a, aspect='auto', norm=LogNorm(vmin=vmin, vmax=vmax), extent=extent)
            ax.text(0.05, 0.95, 'Median: %.4g' % np.median(a), transform=ax.transAxes, va='top', ha='left')

        axs[3, j].set_xlabel(other_col)
        axs[0, j].set_title('gg*_%s' % name)

    for ax in axs[:, 1:].flatten():
        ax.yaxis.set_ticklabels([])

    axs[0, 0].set_ylabel('Mean along %s\n\nTime idx' % (avg_type))
    axs[1, 0].set_ylabel('Std along %s\nof the time diff gains\n\nTime idx' % (avg_type))
    axs[2, 0].set_ylabel('Std along %s\nof the SB diff gains\n\nTime idx' % (avg_type))
    axs[3, 0].set_ylabel('400 s time avg\nstd along %s\nof the SB diff gains\n\nTime idx' % (avg_type))

    cbs = ColorbarSetting(ColorbarInnerPosition(height='70%', pad=0.7))
    cbs.add_colorbar(im_map, axs[0, j])

    fig.tight_layout(pad=0.4)
    fig.savefig(os.path.join(out_dir, '%s_mean_std_%s_gg.pdf' % (prefix, avg_type)))

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

    # plot_sol_mean_and_std(
    #     freqs, gg, pols, args.out_dir, avg_time=True, prefix=args.out_prefix + "sol_DI"
    # )
    # plot_sol_mean_and_std(
    #     freqs, gg, pols, args.out_dir, avg_time=False, prefix=args.out_prefix + "sol_DI"
    # )

    plot_sol_mean_and_std(freqs, gg_stokes, pol_stokes, args.out_dir, avg_type='stations',
                          prefix=args.out_prefix + 'sol_DI', max_station=48)
    plot_sol_mean_and_std(freqs, gg_stokes, pol_stokes, args.out_dir, avg_type='frequency',
                          prefix=args.out_prefix + 'sol_DI')
    # plot_sol_mean_and_std(freqs, gg_stokes, pol_stokes, args.out_dir, avg_type='time',
    #                       prefix=args.out_prefix + 'sol_DI')
    


    # plot_sol_mean_and_std(
    #     freqs,
    #     gg_stokes,
    #     pol_stokes,
    #     args.out_dir,
    #     avg_time=True,
    #     prefix=args.out_prefix + "sol_stokes_DI",
    # )
    # plot_sol_mean_and_std(
    #     freqs,
    #     gg_stokes,
    #     pol_stokes,
    #     args.out_dir,
    #     avg_time=False,
    #     prefix=args.out_prefix + "sol_stokes_DI",
    # )

    print("All done!")


if __name__ == "__main__":
    main()
