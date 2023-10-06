import numpy as np
import matplotlib as mpl

import matplotlib.pyplot as plt

from multiprocessing import Pool

from casacore import tables as tb
from matplotlib import ticker
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable


mpl.style.use("default")
mpl.rcParams["xtick.direction"] = "in"
mpl.rcParams["ytick.direction"] = "in"
mpl.rcParams["xtick.top"] = True
mpl.rcParams["ytick.right"] = True
mpl.rcParams["image.interpolation"] = "none"
mpl.rcParams["image.origin"] = "lower"
mpl.rcParams["image.aspect"] = "auto"


class Flags(object):
    def __init__(self, ms_files):
        self.files = ms_files

    def read_ms(self, filename):
        with tb.table(filename, readonly=True, ack=False) as t:
            flag_col = t.getcol("FLAG")
            time_all = t.getcol("TIME")

            with tb.table(
                filename + "/SPECTRAL_WINDOW", readonly=True, ack=False
            ) as tf:
                freq = tf.getcol("CHAN_FREQ")[0] / 1.0e6  # MHz

            time = np.unique(time_all)
            nbls = int(time_all.size / time.size)

            print(filename, nbls)

            flag_xx = flag_col[..., 0]  # pol XX

            flag_occ = np.zeros((time.size, freq.size))
            for i in range(freq.size):
                row_i = 0
                for j in range(time.size):
                    row_f = row_i + nbls
                    flag_occ[j, i] = np.sum(flag_xx[row_i:row_f, i]) / nbls
                    row_i += nbls

        return (
            flag_occ,
            freq,
            time,
        )

    def plot_flags_occupancy(self, num_cores=9):
        pool = Pool(num_cores)
        data_list = pool.map(self.read_ms, self.files)
        flag_occ_all = np.hstack([d[0] for d in data_list])
        freq_all = np.ravel([d[1] for d in data_list])
        time = np.unique([d[2] for d in data_list])

        # plotting
        flag_occ_all_nozero = np.copy(flag_occ_all)
        flag_occ_all_nozero[flag_occ_all_nozero == 0.0] = 1.0e-15

        formatter = ticker.ScalarFormatter(useMathText=True)
        formatter.set_powerlimits((-3, 3))

        fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))

        max_rfi = np.max(flag_occ_all_nozero * 100) + 1 #+1 to avoid lognorm error when rfi is zero
        im1 = ax1.imshow(
            flag_occ_all_nozero * 100,
            cmap="inferno",  # cmap="RdYlGn_r",
            norm=LogNorm(vmin=.1, vmax=max_rfi),
            extent=[
                freq_all[0],
                freq_all[-1],
                0,
                (time[-1] - time[0]) / 3600,
            ],
        )
        ax1.set_xticklabels(" ")
        ax1.set_ylabel("Obs. Time (hr)")

        axins1 = inset_axes(
            ax1,
            width="100%",
            height="6%",
            loc="upper center",
            bbox_to_anchor=(0, 0.1, 1, 1),
            bbox_transform=ax1.transAxes,
            borderpad=0,
        )
        cb1 = fig.colorbar(
            im1, cax=axins1, orientation="horizontal"
        )  # .set_label("Jy/pixel")
        cb1.ax.xaxis.set_ticks_position("top")
        cb1.ax.xaxis.set_label_position("top")
        cb1.ax.xaxis.set_major_formatter(formatter)
        cb1.set_label("RFI %")

        flag_occ_time = np.sum(flag_occ_all, axis=1) / freq_all.size
        flag_occ_freq = np.sum(flag_occ_all, axis=0) / time.size

        flag_occ_sb = []
        freq_sb = []
        sb_i = 0
        for c in range(int(freq_all.size / 3)):
            sb_f = sb_i + 3
            freq_sb = np.append(freq_sb, np.mean(freq_all[sb_i:sb_f]))
            flag_occ_sb = np.append(flag_occ_sb, np.sum(flag_occ_freq[sb_i:sb_f]) / 3)
            sb_i += 3

        divider = make_axes_locatable(ax1)
        axbottom = divider.append_axes("bottom", size=1.2, pad=0.1)
        axbottom.plot(
            np.sort(freq_sb),
            flag_occ_sb * 100,
            lw=0.8,
            color="black",
            marker=".",
            ms=3,
        )
        axbottom.set_xlabel("Frequency (MHz)")
        axbottom.set_ylabel("RFI %")
        # yrange = axbottom.get_ylim()
        # axbottom.set_ylim(yrange[0], yrange[1])
        axbottom.set_ylim(0, max_rfi)
        axbottom.margins(x=0)
        axbottom.grid(axis="both")

        axright = divider.append_axes("right", size=1.2, pad=0.1)
        axright.plot(
            flag_occ_time * 100,
            (time - time[0]) / 3600,
            lw=0.8,
            color="black",
            marker=".",
            ms=0.01,
        )
        axright.set_yticklabels(" ")
        axright.set_xlabel("RFI %")
        axright.margins(y=0)
        axright.set_xlim(0, max_rfi)
        axright.grid(axis="both")

        plt.tight_layout()
        return fig
