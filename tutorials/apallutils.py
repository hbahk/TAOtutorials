# -*- coding: utf-8 -*-
"""
@Author: Hyeonguk Bahk
@Email: bahkhyeonguk@gmail.com
@Date: 2024-05-21

@Filename: apallutils.py
@Brief: This file contains the utility functions for the aperture extraction.
@Details: This includes the functions for the aperture tracing and summing in
the context of the astronomical spectroscopy.
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.stats import sigma_clip, gaussian_fwhm_to_sigma
from astropy.visualization import ZScaleInterval
from astropy import units as u
from ccdproc import cosmicray_lacosmic, gain_correct
from numpy.polynomial.chebyshev import chebfit, chebval
from scipy.ndimage import median_filter
from scipy.signal import peak_widths
from skimage.feature import peak_local_max
from photutils.centroids import centroid_com


def fit_background(x, apall, fitmask, sigma_sigclip, order_skymodel):
    x_sky = x[fitmask]
    apall_sky = apall[fitmask]
    # sigma clipping
    mask_sigclip = sigma_clip(apall_sky, sigma=sigma_sigclip, masked=True).mask

    # fitting the sky background
    coeffs_sky, fitfull_sky = chebfit(
        x_sky[~mask_sigclip],
        apall_sky[~mask_sigclip],
        deg=order_skymodel,
        full=True,
    )

    # sky background model
    x = np.arange(apall.size)
    skymodel = chebval(x, coeffs_sky)

    # calculating the rms
    residual_skymodel = fitfull_sky[0][0]
    rms_skymodel = np.sqrt(residual_skymodel / len(x_sky[~mask_sigclip]))

    return skymodel, rms_skymodel, mask_sigclip


def get_skymask(apall, sky_limit_peak, center=None):
    if center is None:
        peak_pix = peak_local_max(
            apall,
            num_peaks=1,
            min_distance=10,
            threshold_abs=np.mean(apall),
        )[0][0]
    else:
        peak_pix = round(center)

    mask_sky = np.ones_like(apall, dtype=bool)
    mask_sky[peak_pix - sky_limit_peak : peak_pix + sky_limit_peak] = False

    return mask_sky, peak_pix


def find_center(data, lcut, rcut, sky_limit_peak, sigma_sigclip, order_skymodel):
    apall = np.sum(data[lcut:rcut, :], axis=0)
    x = np.arange(apall.size)

    # mask for the sky region and initial peak pixel
    mask_sky, peak_pix = get_skymask(apall, sky_limit_peak)

    # fitting the sky background
    skymodel, rms_skymodel, mask_sigclip = fit_background(
        x,
        apall,
        mask_sky,
        sigma_sigclip=sigma_sigclip,
        order_skymodel=order_skymodel,
    )

    apall_skysub = apall - skymodel

    # center of the source
    center_com = centroid_com(apall_skysub)[0]

    # fwhm of the source
    pwresult = peak_widths(apall_skysub, [round(peak_pix)], rel_height=0.5)
    fwhm = pwresult[0][0]

    return center_com, fwhm


def aperture_sum(
    cut, center, apsum_sigma, ap_sigma, sky_limit_peak, sigma_sigclip, order_skymodel, rdnoise
):
    x = np.arange(len(cut))

    mask_sky, peak_pix = get_skymask(cut, sky_limit_peak, center)
    # aperture size = trace center +/- apsum_sigma * ap_sigma
    appos_lower = int(center - apsum_sigma * ap_sigma)
    appos_upper = int(center + apsum_sigma * ap_sigma)

    # fitting the sky background
    skymodel, rms_skymodel, mask_sigclip = fit_background(
        x,
        cut,
        mask_sky,
        sigma_sigclip=sigma_sigclip,
        order_skymodel=order_skymodel,
    )

    # subtracting the sky background
    cut_skysub = cut - skymodel
    source = cut_skysub[appos_lower:appos_upper]

    # aperture sum
    ap_sum = np.sum(source)

    # calculating the uncertainty
    sig = rdnoise**2 + source + skymodel[appos_lower:appos_upper]
    ap_sig = np.sqrt(np.sum(sig))

    return ap_sum, ap_sig


class InstrumentalSpectrum:
    def __init__(
        self,
        fpath,
        sky_limit_peak,
        sigma_sigclip,
        order_skymodel,
        order_aptrace,
        sigma_clip_aptrace,
        slice_width_aptrace,
        apsum_sigma,
        trim=[0, 4220, 670, 810],
    ):
        self.fname = fpath
        self.TRIM = trim
        self.SKY_LIMIT_PEAK = sky_limit_peak
        self.SIGMA_SIGCLIP = sigma_sigclip
        self.ORDER_SKYMODEL = order_skymodel
        self.ORDER_APTRACE = order_aptrace
        self.SIGMA_CLIP_APTRACE = sigma_clip_aptrace
        self.SLICE_WIDTH_APTRACE = slice_width_aptrace
        self.APSUM_SIGMA = apsum_sigma

        self.hdu = fits.open(fpath)
        self.data = self.hdu[0].data[
            self.TRIM[0] : self.TRIM[1], self.TRIM[2] : self.TRIM[3]
        ]
        self.hdr = self.hdu[0].header
        self.EXPTIME = self.hdr["EXPTIME"]
        self.RDNOISE = self.hdr["RDNOISE"]

        self.ccd = CCDData(data=self.data, header=self.hdr, unit="adu")
        self.gcorr = gain_correct(self.ccd, gain=self.hdr["GAIN"] * u.electron / u.adu)
        self.crrej = cosmicray_lacosmic(
            self.gcorr,
            sigclip=7,
            readnoise=self.hdr["RDNOISE"] * u.electron,
            verbose=False,
        )

    def set_aptrace(self):
        slice_width = self.SLICE_WIDTH_APTRACE
        lcut, rcut = 0, slice_width
        midpoints = []
        aptrace = []
        aptrace_fwhm = []
        while rcut < len(self.crrej.data):
            rcut = lcut + slice_width
            midpoints.append((lcut + rcut) / 2)
            center, fwhm = find_center(
                self.crrej.data,
                lcut,
                rcut,
                self.SKY_LIMIT_PEAK,
                self.SIGMA_SIGCLIP,
                self.ORDER_SKYMODEL,
            )
            aptrace.append(center)
            aptrace_fwhm.append(fwhm)
            lcut += slice_width
            rcut += slice_width

        self.aptrace = np.array(aptrace)
        self.aptrace_fwhm = np.array(aptrace_fwhm)
        self.aptrace_sigma = self.aptrace_fwhm * gaussian_fwhm_to_sigma
        midpoints = np.array(midpoints)
        self.aptrace_midpoints = midpoints

        # fitting the trace line
        coeffs_aptrace_init = chebfit(midpoints, self.aptrace, deg=self.ORDER_APTRACE)

        # sigma clipping
        mask_aptrace = sigma_clip(
            self.aptrace - chebval(midpoints, coeffs_aptrace_init),
            sigma=self.SIGMA_CLIP_APTRACE,
            masked=True,
        ).mask
        self.mask_aptrace = mask_aptrace

        # fitting the trace line again
        self.coeffs_aptrace = chebfit(
            midpoints[~mask_aptrace],
            self.aptrace[~mask_aptrace],
            deg=self.ORDER_APTRACE,
        )
        self.residuals_aptrace = self.aptrace - chebval(midpoints, self.coeffs_aptrace)
        self.rejected_aptrace = ~np.in1d(midpoints, midpoints[~mask_aptrace])
        self.rms_aptrace = np.sqrt(np.mean(self.residuals_aptrace[~mask_aptrace] ** 2))

        # aperture trace model
        ximg = np.arange(self.crrej.data.shape[0])
        self.aptrace_model = chebval(ximg, self.coeffs_aptrace)
        self.aptrace_model_midpoints = chebval(midpoints, self.coeffs_aptrace)

    def set_apsum(self):
        ap_sum = []
        ap_sig = []
        self.ap_fwhm = np.median(self.aptrace_fwhm[~self.mask_aptrace])
        self.ap_sigma = self.ap_fwhm * gaussian_fwhm_to_sigma

        for i in range(self.crrej.data.shape[0]):
            cut_i = self.crrej.data[i, :]
            center_i = self.aptrace_model[i]

            ap_sum_i, ap_sig_i = aperture_sum(
                cut_i,
                center_i,
                self.APSUM_SIGMA,
                np.median(self.ap_sigma),
                self.SKY_LIMIT_PEAK,
                self.SIGMA_SIGCLIP,
                self.ORDER_SKYMODEL,
                self.RDNOISE,
            )
            ap_sum.append(ap_sum_i)
            ap_sig.append(ap_sig_i)

        self.ap_sum = np.array(ap_sum) / self.EXPTIME
        self.ap_sig = np.array(ap_sig) / self.EXPTIME
        
    def draw_img(self):
        fig = plt.figure(figsize=(12, 6))
        ax = fig.add_subplot(111)
        interval = ZScaleInterval()
        vmin, vmax = interval.get_limits(self.crrej.data)
        im = ax.imshow(self.crrej.data.T, cmap="gray", origin="lower", vmin=vmin, vmax=vmax)
        ax.set_title("Cosmic Ray Removed Feige34 Image")
        fig.show()
        return fig

    def plot_aptrace(self):
        fig, axes = plt.subplots(
            2, 1, figsize=(12, 6), gridspec_kw={"height_ratios": [3, 1]}
        )
        ax = axes[0]
        interval = ZScaleInterval()
        vmin, vmax = interval.get_limits(self.crrej.data)
        ax.imshow(
            self.crrej.data.T,
            cmap="gray",
            origin="lower",
            vmin=vmin,
            vmax=vmax,
            aspect="auto",
        )
        ax.scatter(
            self.aptrace_midpoints,
            self.aptrace,
            marker="+",
            c="dodgerblue",
            s=100,
            lw=2,
            label="Data",
        )
        ax.scatter(
            self.aptrace_midpoints[self.rejected_aptrace],
            self.aptrace[self.rejected_aptrace],
            marker="x",
            c="salmon",
            s=100,
            lw=1,
            label="Rejected",
        )
        ximg = np.arange(self.crrej.data.shape[0])
        ax.plot(
            ximg, self.aptrace_model, c="tab:red", lw=1, label="Chebyshev Trace Model"
        )
        ax.plot(
            self.aptrace_midpoints,
            self.aptrace_model_midpoints + 5 * self.aptrace_sigma,
            c="r",
            lw=0.8,
            ls="--",
            label="5-sigma",
        )
        ax.plot(
            self.aptrace_midpoints,
            self.aptrace_model_midpoints - 5 * self.aptrace_sigma,
            c="r",
            lw=0.8,
            ls="--",
        )
        ax.plot(
            self.aptrace_midpoints,
            self.aptrace_model_midpoints + 10 * self.aptrace_sigma,
            c="r",
            lw=0.8,
            ls=":",
            label="10-sigma",
        )
        ax.plot(
            self.aptrace_midpoints,
            self.aptrace_model_midpoints - 10 * self.aptrace_sigma,
            c="r",
            lw=0.8,
            ls=":",
        )
        ax.plot(
            self.aptrace_midpoints,
            self.aptrace_model_midpoints + self.APSUM_SIGMA * np.median(self.aptrace_sigma),
            c="g",
            lw=0.8,
            ls="-",
            label="Aperture Sum Region",
        )
        ax.plot(
            self.aptrace_midpoints,
            self.aptrace_model_midpoints - self.APSUM_SIGMA * np.median(self.aptrace_sigma),
            c="g",
            lw=0.8,
            ls="-",
        )
        ax.set_title("Aperture Trace")
        ax.set_ylabel("Pixel Number (Spatial Direction)")
        ax.set_xticks([])
        ax.legend()

        resax = axes[1]
        resax.plot([0, ximg[-1]], [0, 0], c="k", ls="--", lw=0.8)
        resax.scatter(
            self.aptrace_midpoints,
            self.residuals_aptrace,
            marker="+",
            c="dodgerblue",
            s=100,
            lw=2,
        )
        resax.scatter(
            self.aptrace_midpoints[self.rejected_aptrace],
            self.residuals_aptrace[self.rejected_aptrace],
            marker="x",
            c="salmon",
            s=100,
            lw=1,
        )
        resax.set_xlabel("Pixel Number (Dispersion Direction)")
        resax.set_ylabel("Residuals [pix]")
        resax.tick_params(axis="x", direction="in", top="on")
        resax_ylim = np.max(np.abs(self.residuals_aptrace)) * np.array([-1, 1])
        resax.set_ylim(resax_ylim)
        resax.set_xlim(0, ximg[-1])
        resax.text(
            0.99,
            0.95,
            f"RMS={self.rms_aptrace:.2f}",
            transform=resax.transAxes,
            ha="right",
            va="top",
        )

        fig.subplots_adjust(hspace=0)
        return fig