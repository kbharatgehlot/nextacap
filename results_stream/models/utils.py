# Random util functions
import numpy as np
import astropy.units as u
import matplotlib as mpl
import matplotlib.pyplot as plt
import streamlit as st
from astropy.io import fits
from astropy.wcs import WCS
from pdf2image import convert_from_path


mpl.style.use("default")
mpl.rcParams["xtick.direction"] = "in"
mpl.rcParams["ytick.direction"] = "in"
mpl.rcParams["xtick.top"] = True
mpl.rcParams["ytick.right"] = True
mpl.rcParams["image.interpolation"] = "none"
mpl.rcParams["image.origin"] = "lower"
mpl.rcParams["image.aspect"] = "auto"

################################################################################################
# All sky maps plotting
###############################################################################################


@st.cache_data
def read_fits_image(image_file, mask_horizon=True):
    hdu = fits.open(image_file)[0]
    data = hdu.data[0, 0, :, :]
    header = hdu.header
    wcs = WCS(header)
    wcs = wcs.dropaxis(2)
    wcs = wcs.dropaxis(2)

    if mask_horizon:
        nx, ny = np.meshgrid(
            np.arange(-0.5 * header["NAXIS1"], 0.5 * header["NAXIS1"])
            * np.abs(header["CDELT1"]),
            np.arange(-0.5 * header["NAXIS2"], 0.5 * header["NAXIS2"])
            * np.abs(header["CDELT2"]),
        )

        horizon_mask = np.zeros_like(data, dtype=bool)
        mask_idx = np.sqrt(np.radians(nx) ** 2 + np.radians(ny) ** 2) >= 1.0
        horizon_mask[mask_idx] = True

        data = np.ma.array(data, mask=horizon_mask)
    return data, wcs


# @st.cache_data(hash_funcs={matplotlib.figure.Figure: lambda _: None})
def all_sky_map(image_file):
    fig = plt.figure()  # figsize=(12,10)
    image_data, wcs = read_fits_image(image_file)
    ax = fig.add_axes([0.1, 0.08, 0.85, 0.85], projection=wcs)
    ax.coords[0].set_auto_axislabel(False)
    ax.coords[1].set_auto_axislabel(False)
    ax.coords[0].set_ticks([0, 90, 180, 270] * u.deg)
    ax.coords[0].set_major_formatter("hh")
    im = ax.imshow(image_data, cmap="Spectral_r", vmin=-0.005, vmax=0.005)
    ax.set_xlabel("RA", fontsize=15)
    ax.set_ylabel("Dec", fontsize=15)
    cbar = fig.colorbar(im, label="Jy/Bm")
    cbar.set_label(label="Jy/Bm", size=15)
    # ax.set_title(f'{obsid}', fontsize=18)
    ax.coords.grid(color="black", alpha=0.7, ls="dotted")
    fig.patch.set_facecolor("white")

    return fig


@st.cache_data
def convert_pdf_to_jpeg(pdf_name):
    image = convert_from_path(
        pdf_name,
        dpi=500,
        fmt="jpeg",
    )
    return image
