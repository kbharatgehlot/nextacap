# Random util functions

import astropy.units as u
import matplotlib.pyplot as plt
import streamlit as st
from astropy.io import fits
from astropy.wcs import WCS
from pdf2image import convert_from_path


################################################################################################
# All sky maps plotting
###############################################################################################


@st.cache_data
def read_fits_image(image_file):
    hdu = fits.open(image_file)[0]
    image_data = hdu.data[0, 0, :, :]
    wcs = WCS(hdu.header)
    wcs = wcs.dropaxis(2)
    wcs = wcs.dropaxis(2)
    return image_data, wcs


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
