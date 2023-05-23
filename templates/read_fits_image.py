"""
Thanks to Nunhokee Ridhima
"""
from astropy.io import fits
from astropy import wcs
from astropy.visualization.wcsaxes import WCSAxes
from astropy.modeling import models, fitting
import numpy as np
import warnings
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys

from argparse import ArgumentParser


parser = ArgumentParser("read and plot fits images")
parser.add_argument(
    "-i",
    "--fits_image",
    help="Fits image path",
    dest="fits_image",
    default="",
)


def main(argv):
    args = parser.parse_args(argv)
    img_obj = Image(args.fits_image)
    img_obj.plot_image(
        vmin=-1e-2,
        vmax=1e-2,
        select=True,
        srcpos=(0, 90),
        beam_const=336,
        cmap="Spectral_r",
        save=True,
    )
    return


class Image(object):
    def __init__(self, fitspath, pix_box=[100, 100]):
        self.fitspath = fitspath
        self.pix_box = pix_box

        with fits.open(self.fitspath) as hdus:
            img_hdu = hdus["PRIMARY"]

            self.data_array = img_hdu.data.squeeze()
            self.header = img_hdu.header
            try:
                self.image_ID = img_hdu.header["OBJECT"]
            except KeyError:
                image_id = int("".join(filter(str.isdigit, fitspath)))
                self.image_ID = image_id
            self.obsdate = img_hdu.header["DATE-OBS"]
            self.image_size = [img_hdu.header["NAXIS1"], img_hdu.header["NAXIS2"]]
            self.xcellsize = np.abs(img_hdu.header["CDELT1"])
            self.ycellsize = np.abs(img_hdu.header["CDELT2"])
            self.beam_major = img_hdu.header["BMAJ"]
            self.beam_minor = img_hdu.header["BMIN"]
            self.beam_parallactic_angle = img_hdu.header["BPA"]
            self.beam_major_px = self.beam_major / self.xcellsize
            self.beam_minor_px = self.beam_minor / self.ycellsize
            self.beam_area = self.beam_major * self.beam_minor * np.pi
            self.beam_npix = self.beam_area / (self.xcellsize * self.ycellsize)
            self.beam_radius_px = np.sqrt(
                self.beam_major_px**2 + self.beam_minor_px**2
            )
            self.mean = np.nanmean(self.data_array)
            self.rms = np.sqrt(np.nanmean(self.data_array**2))
            self.std = np.nanstd(self.data_array)
            self.polarization = img_hdu.header["CRVAL4"]
            region = self.data_array[0 : self.pix_box[0], 0 : self.pix_box[1]]
            self.mean_across_box = np.nanmean(region)
            self.std_across_box = np.nanstd(region)
            self.rms_across_box = np.sqrt(np.nanmean(region**2))

    def src_pix(self, src_pos):
        ra, dec = src_pos
        w = wcs.WCS(self.fitspath)
        xpix, ypix = w.all_world2pix(ra, dec, 1, 0, 0)[:2]
        if (xpix > 0) and (xpix < self.image_size[0]):
            if (ypix > 0) and (ypix < self.image_size[0]):
                xpix, ypix = np.around(xpix, 0).astype(int), np.around(ypix, 0).astype(
                    int
                )
        else:
            xpix, ypix = np.nan, np.nan

        return xpix, ypix

    def src_flux(self, src_pos, beam_const=1, deconvol=False):
        if deconvol:
            pflux, tflux, std = self.fit_gaussian(src_pos, beam_const)
        else:
            select = self._select_region(src_pos, beam_const)
            if select is None:
                pflux, tflux, std = np.nan, np.nan, np.nan
            else:
                region = self.data_array.squeeze()[select]
                pflux = np.nanmax(np.abs(region))
                # might have negative pixels as well
                npoints = len(np.where(select == True)[0])
                tflux = np.nansum(np.abs(region)) / npoints
                std = np.nanstd(region)

        return pflux, tflux, std

    def _select_region(self, src_pos, beam_const):
        ra_pix, dec_pix = self.src_pix(src_pos)
        if np.isnan(ra_pix) or np.isnan(dec_pix):
            warnings.warn("Source position out of bounds")
            select = None
        else:
            l_axis = np.arange(self.image_size[0])
            m_axis = np.arange(self.image_size[1])
            ll, mm = np.meshgrid(l_axis, m_axis)
            R = np.sqrt((ll - ra_pix) ** 2 + (mm - dec_pix) ** 2)
            select = R < beam_const * self.beam_radius_px

        return select

    def fit_gaussian(self, src_pos, bm_const):
        select = self._select_region(src_pos, bm_const)
        if select is None:
            gauss_pflux, gauss_tflux, gauss_std = np.nan, np.nan, np.nan
        else:
            l_axis = np.arange(self.image_size[0])
            m_axis = np.arange(self.image_size[1])
            ll, mm = np.meshgrid(l_axis, m_axis)
            gauss_data = self.data_array
            gauss_data[~select] = 0
            gauss_data = gauss_data.reshape(self.image_size)
            peak_val = np.nanmax(gauss_data)
            inds = np.where(gauss_data == peak_val)
            mod = models.Gaussian2D(
                peak_val,
                inds[1],
                inds[0],
                self.beam_major_px / 2,
                self.beam_minor_px / 2,
                theta=self.beam_parallactic_angle * np.pi / 180,
            )
            fit_p = fitting.LevMarLSQFitter()
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
            try:
                mod = models.Gaussian2D(
                    peak_val,
                    inds[1],
                    inds[0],
                    self.beam_major_px / 2,
                    self.beam_minor_px / 2,
                    theta=self.beam_parallactic_angle * np.pi / 180,
                )
                gauss_mod = fit_p(mod, ll, mm, gauss_data)
            except ValueError:
                mod = models.Gaussian2D(
                    peak_val,
                    inds[0][1],
                    inds[0][0],
                    self.beam_major_px / 2,
                    self.beam_minor_px / 2,
                    theta=self.parallactic_angle * np.pi / 180,
                )
                gauss_mod = fit_p(mod, ll, mm, gauss_data)
            gauss_pflux = gauss_mod.amplitude.value
            npoints = len(np.where(select == True)[0])
            gauss_tflux = np.nansum(gauss_mod(ll, mm)) / npoints
            gauss_std = np.nanstd(gauss_mod(ll, mm))
        return gauss_pflux, gauss_tflux, gauss_std

    def plot_image(
        self,
        select=False,
        srcpos=None,
        beam_const=None,
        vmin=None,
        vmax=None,
        cmap="viridis",
        save=None,
        figname=None,
        dpi=200,
    ):
        my_wcs = wcs.WCS(self.header, naxis=[wcs.WCSSUB_CELESTIAL])
        fig = plt.figure(figsize=(12, 8), dpi=dpi)
        ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=my_wcs)
        ax2 = fig.add_axes(ax)
        if vmin is None:
            vmax = np.nanmin(self.data_array)
        if vmax is None:
            vmax = np.nanmax(self.data_array)
        if select:
            plot_data = self.data_array.squeeze() * self._select_region(
                srcpos, beam_const
            )
        else:
            plot_data = self.data_array.squeeze()
        im = ax2.imshow(plot_data, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)

        # change colobar labels to mJy
        func = lambda x, pos: "{:g}".format(x * 1000)
        fmt = mpl.ticker.FuncFormatter(func)

        cbar = plt.colorbar(im, ax=ax, format=fmt, extend="both")
        cbar.set_label("mJy/PSF")  # ( f"m{self.header['BUNIT'].title()}")
        ax2.set_xlabel("R.A. [deg]")
        ax2.set_ylabel("Dec [deg]")
        plt.grid(ls="dotted", color="black")

        # TODO: Fix these hacks
        cx, cy = plot_data.shape[0] // 2, plot_data.shape[1] // 2
        circle = plt.Circle((cx, cy), 450, color="k", fill=False, lw=2)
        ax.add_patch(circle)

        casAx, casAy = plot_data.shape[0] // 2 + 150, plot_data.shape[1] // 4
        circle = plt.Circle((casAx, casAy), 250, color="k", fill=False, lw=2)
        ax.add_patch(circle)

        cx, cy = cx - 60, cy - 60
        rect = plt.Rectangle((cx, cy), 120, 120, color="k", fc="none", lw=2)
        ax.add_patch(rect)

        fig.tight_layout()
        if save:
            if figname is None:
                figname = self.fitspath.replace(".fits", ".png")
            plt.savefig(figname, dpi=dpi, bbox_inches="tight")
            plt.close()
        else:
            plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])
