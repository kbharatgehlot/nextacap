import sys
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter

import matplotlib.pyplot as plt
from pspipe import settings, database

import matplotlib as mpl

mpl.rcParams["image.cmap"] = "Spectral_r"
mpl.rcParams["image.origin"] = "lower"
mpl.rcParams["image.interpolation"] = "nearest"
mpl.rcParams["axes.grid"] = True

parser = ArgumentParser(
    description="Plot DD calibration solutions", formatter_class=RawTextHelpFormatter
)

parser.add_argument("toml", help="config file")
parser.add_argument("--obsid", help="obsid")
parser.add_argument("--outdir", help="outdir")

args = parser.parse_args(sys.argv[1:])


s = settings.Settings.load_with_defaults(args.toml)
rev = database.VisRevision(s)

data = rev.get_data(args.obsid)
# data.do_flag()


ps_gen = data.get_ps_gen(
    rmean_freqs=True, umax=250, umin=50, du=8, window_fct="hann", eor_bin_name="2"
)
# ps_gen.get_ps2d(data.i).plot(vmin=2e1, vmax=1e4)

fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(9, 5), sharey=True, dpi=200)
ps_gen.get_ps2d(data.i).plot(
    ax=ax1, title="PS2D Stokes I", vmin=2e1, vmax=1e6, wedge_lines=[90], z=ps_gen.eor.z
)
ps_gen.get_ps2d(data.v).plot(
    ax=ax2, title="PS2D Stokes V", vmin=2e1, vmax=1e6, wedge_lines=[90], z=ps_gen.eor.z
)
fig.tight_layout()
plt.savefig(f"{args.outdir}/{args.obsid}_IV_power_spectrum.png", bbox_inches="tight")
plt.close()
