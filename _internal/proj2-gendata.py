#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 19 16:40:23 2023

Data makeup for TAO tutorial - PSF photometry

@author: Hyeonguk Bahk
"""

from pathlib import Path

from astropy import units as u
from astropy.io import fits
from astropy.visualization import ZScaleInterval
from astropy.wcs import WCS
from matplotlib import pyplot as plt

from exgalcosutils.from_legacy_survey import dr2_rgb, dr10_griz_rgb

HOME = Path.home()
DATADIR = HOME/'class'/'tao23'/'TAO23'/'tutorials'/'data'/'proj2'
# plt.rcParams["figure.figsize"] = (7, 5)
plt.rcParams["font.family"] = 'Times New Roman'
plt.rcParams["font.size"] = 15
plt.rcParams["text.usetex"] = False
plt.rcParams["mathtext.fontset"] = 'cm'
#%%
# cra, cdec = 92.8804, -69.1214
cra, cdec = 92.8290, -69.1046

# size = 2*u.arcmin
size = 1*u.arcmin
pix_scale = 0.262

if type(size) == u.quantity.Quantity:
    SCALE = 25/9*1e-4*u.deg * pix_scale  # scale angle of a pixel
    pix_size = round(((2*size/SCALE).decompose().value))
elif type(size) == int:
    pix_size = size
else:
    raise ValueError('`size` should be either type of int '+
                     'astropy Quantity in angular unit')

img_query_url = 'https://www.legacysurvey.org/viewer/cutout.fits?'\
+f'ra={cra:.8f}&dec={cdec:.8f}&width={pix_size}&height={pix_size}'\
    + f'&layer=ls-dr10&pixscale={pix_scale}&bands=griz'

img_hdu = fits.open(img_query_url)
h = img_hdu[0].header
wcs = WCS(h)

if h['BANDS'] == 'grz':
    gimg = img_hdu[0].data[0]
    rimg = img_hdu[0].data[1]
    zimg = img_hdu[0].data[2]
    rgbimg = dr2_rgb(rimgs=[gimg, rimg, zimg], bands=['g','r','z'])
elif h['BANDS'] == 'griz': 
    gimg = img_hdu[0].data[0]
    rimg = img_hdu[0].data[1]
    iimg = img_hdu[0].data[2]
    zimg = img_hdu[0].data[3]
    rgbimg = dr10_griz_rgb(rimgs=[gimg, rimg, iimg, zimg],
                           bands=['g','r','i','z'])

img_hdu.writeto(DATADIR/'NGC2210.grz.fits', overwrite=True)
#%%
cmap = 'gray'
minpad = 0.3

interval = ZScaleInterval()
gn, gx = interval.get_limits(gimg)
rn, rx = interval.get_limits(rimg)
zn, zx = interval.get_limits(zimg)


fig = plt.figure(figsize=(20,5))
ax0 = fig.add_subplot(141, projection=wcs, slices=('x','y',0))
ax1 = fig.add_subplot(142, projection=wcs, slices=('x','y',0))
ax2 = fig.add_subplot(143, projection=wcs, slices=('x','y',0))
ax3 = fig.add_subplot(144, projection=wcs, slices=('x','y',0))
ax0.imshow(rgbimg)
ax1.imshow(gimg, vmin=gn, vmax=gx, cmap=cmap)
ax2.imshow(rimg, vmin=rn, vmax=rx, cmap=cmap)
ax3.imshow(zimg, vmin=zn, vmax=zx, cmap=cmap)

fig.subplots_adjust(wspace=0.0)

lon = ax0.coords['ra']
lat = ax0.coords['dec']
lon.set_axislabel('RA', minpad=minpad)
lat.set_axislabel('Dec')

bbox = dict(boxstyle='round', facecolor='w', alpha=0.9)

for ax, band in zip([ax1, ax2, ax3], ['g','r','z']):
    lon = ax.coords['ra']
    lat = ax.coords['dec']
    lat.set_ticks_visible(False)
    lat.set_ticklabel_visible(False)
    lon.set_axislabel('RA', minpad=minpad)
    
    ax.text(0.98, 0.02, f'${band}$-band', transform=ax.transAxes,
            va='bottom', ha='right', c='k', bbox=bbox)

ax0.text(0.98, 0.02, '$grz$ false color', transform=ax0.transAxes,
         va='bottom', ha='right', c='k', bbox=bbox)
fig.suptitle('NGC 2210 - DESI Legacy Survey $grz$ image')

#%%

import numpy as np
from time import time, ctime
import warnings
from astropy.table import Table
from astropy.modeling.fitting import LevMarLSQFitter, LMLSQFitter
from astropy.stats import sigma_clipped_stats, gaussian_fwhm_to_sigma
from astropy.table import Table, vstack, hstack
from astropy.nddata import NDData, CCDData, Cutout2D
from astropy.utils.exceptions import AstropyWarning
from astropy.visualization import simple_norm
from photutils.aperture import CircularAperture, ApertureStats
from photutils.detection import DAOStarFinder, find_peaks
from photutils.psf import extract_stars, EPSFBuilder
from photutils.psf import IterativelySubtractedPSFPhotometry
from photutils.psf import DAOPhotPSFPhotometry, BasicPSFPhotometry
from photutils.detection import IRAFStarFinder
from photutils.psf import DAOGroup, DBSCANGroup, IntegratedGaussianPRF
from photutils.background import MMMBackground, MADStdBackgroundRMS
from photutils.segmentation import make_source_mask
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import patches
import sep

warnings.simplefilter('ignore', category=AstropyWarning)

HOME = Path.home()
DATADIR = HOME/'class'/'tao23'/'TAO23'/'tutorials'/'data'/'proj2'
#%%
def zscale_imshow(ax, img, vmin=None, vmax=None, **kwargs):
    if vmin==None or vmax==None:
        interval = ZScaleInterval()
        _vmin, _vmax = interval.get_limits(img)
        if vmin==None:
            vmin = _vmin
        if vmax==None:
            vmax = _vmax
    im = ax.imshow(img, vmin=vmin, vmax=vmax, origin='lower', **kwargs)
    
    return im

def read_sci_data(fpath, show=True, newbyteorder=False, bkgsubtract=False):
    ccd_raw = CCDData.read(fpath, unit='nmgy') # be careful about the unit!
    ccd = CCDData.read(fpath, unit='nmgy')
    
    data = ccd.data
    
    data_mean, data_med, data_std = sigma_clipped_stats(data, sigma=2.)
    ccd.mask = data < data_mean - 3*data_std 
    ccd.data[ccd.mask] = data_med
    ccd.data -= data_med
    ccd.mask = (ccd.mask) & (data > 10.)
    if newbyteorder:
        ccd.data = ccd.data.byteswap().newbyteorder()
        
    n = len(ccd.data)
    if bkgsubtract:
        bkgsubs = []
        for i in range(n):
            bkgsub = background_substraction(ccd.data[i])
            bkgsubs.append(bkgsub)
        ccd.data = np.array(bkgsubs)
    
    if show:
        fig, axes = plt.subplots(1, n, figsize=(5*n, 5))
        for i in range(n):
            im = zscale_imshow(axes[i], ccd.data[i], cmap='gray')
            divider = make_axes_locatable(axes[i])
            cax = divider.append_axes('right', size='5%', pad=0.05)
            fig.colorbar(im, cax=cax, orientation='vertical')
            
            band = ccd.header[f'BAND{i}']
            
            axes[i].text(0.98, 0.02, f'${band}$-band',
                         transform=axes[i].transAxes,
                         va='bottom', ha='right', c='k',
                         bbox=dict(boxstyle='round', facecolor='w', alpha=0.9))
            
        plt.tight_layout()
        
    return ccd, ccd_raw

def background_substraction(img, box=100, filt=3, show=True):
    mask = make_source_mask(img, nsigma=2, npixels=1, dilate_size=3)
    bkg_sep = sep.Background(img.byteswap().newbyteorder(), 
                             mask=mask, bw=box, bh=box, fw=filt, fh=filt)
    bkgsub_sep = img - bkg_sep.back()
    
    if show:
        fig, axs = plt.subplots(2, 3, figsize=(10, 8))
        axs[1, 1].axis("off")
        
        data2plot = [
            dict(ax=axs[0, 0], arr=img,              title="Original data"),
            dict(ax=axs[0, 1], arr=bkg_sep.back(),   title=f"bkg (filt={filt:d}, box={box:d})"),
            dict(ax=axs[0, 2], arr=bkgsub_sep,       title="bkg subtracted"),
            dict(ax=axs[1, 0], arr=mask,             title="Mask"),
            # dict(ax=axs[1, 1], arr=bkg_sep.background_mesh,   title="bkg mesh"),
            dict(ax=axs[1, 2], arr=10*bkg_sep.rms(), title="10 * bkg RMS")
        ]
        
        for dd in data2plot:
            im = zscale_imshow(dd['ax'], dd['arr'])
            cb = fig.colorbar(im, ax=dd['ax'], orientation='horizontal')
            cb.ax.set_xticklabels(cb.get_ticks().astype(int), rotation=45)
            dd['ax'].set_title(dd['title'])
        
        plt.tight_layout()
    
    return bkgsub_sep
#%%

ccd, ccd_raw = read_sci_data(DATADIR/'NGC2210.grz.fits',bkgsubtract=True)

#%%
bands = ['g', 'r', 'z']
i = 0
band = bands[i]
data = ccd.data[i]
mask = ccd.mask[i]
find_mask = np.zeros_like(data, dtype=bool)
# find_mask[250:650, 250:650] = True
ll, hh = 0, 200
find_mask[ll:hh, ll:hh] = True
peaks_tbl = find_peaks(data, threshold=0.05, mask=find_mask)
peaks_tbl['peak_value'].info.format = '%.8g'

print(peaks_tbl)

# select stars within the cutout of specified size
size = 25
hsize = (size - 1) / 2
x = peaks_tbl['x_peak']  
y = peaks_tbl['y_peak']  
bound = ((x > hsize) & (x < (data.shape[1] -1 - hsize)) &
         (y > hsize) & (y < (data.shape[0] -1 - hsize)))
bright = (peaks_tbl['peak_value'] > 0.5) & (peaks_tbl['peak_value'] < 10.)

bbmask = (bound & bright)
bx, by = x[bbmask], y[bbmask]

isolated = [False if np.count_nonzero(np.sqrt((bx-xi)**2+(by-yi)**2)<size) > 1
            else True for xi, yi in zip(bx, by)]

mask_stars = isolated
mask_stars = np.ones_like(bx, dtype=bool)

stars_tbl = Table()
stars_tbl['x'] = bx[mask_stars]  
stars_tbl['y'] = by[mask_stars]

print(stars_tbl)

nddata = NDData(data=data, mask=ccd.mask[i])

# plot the locations of our selected stars
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
im = zscale_imshow(ax, data, cmap='gray')
ax.plot(stars_tbl['x'], stars_tbl['y'], '.r')
rect = patches.Rectangle((ll,ll), hh-ll, hh-ll, ec='tab:red', ls='--',
                         hatch='/', fc='none', lw=1)
ax.add_patch(rect)

#%%
# extract cutouts of our selected stars
stars = extract_stars(nddata, stars_tbl, size=size)

# plot cutouts
nrows = 5
ncols = 5
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, 10),
                       squeeze=True)
ax = ax.ravel()
for i in range(nrows * ncols):
    norm = simple_norm(stars[i], 'log', percent=99.0)
    ax[i].imshow(stars[i], norm=norm, origin='lower', cmap='viridis')

# derive fwhm from aperture statistics
cap = CircularAperture(stars.center_flat, size/5.)
apstat = ApertureStats(data, cap)

fwhm_mean, fwhm_med, fwhm_std = sigma_clipped_stats(apstat.fwhm.value,
                                                    sigma=3.0)
#%%

def get_epsf(stars, band,
             oversampling=4, maxiters=10, smoothing_kernel='quadratic',
             show=True, ncols=10, figsize=(10, 6)):
    
    if show:
        nrows = int(np.ceil(len(stars)/ncols))
        fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize,
                               squeeze=True)
        axs = axs.ravel()
        for ax in axs:
            ax.axis('off')
        for i in range(len(stars)):
            norm = simple_norm(stars[i], 'log', percent=99.)
            axs[i].imshow(stars[i], norm=norm, origin='lower', cmap='viridis')
        fig.text(0.95,0.05,f'{band}-band',
                  ha='right', va='bottom', fontsize=15)
        # plt.title(f'Picked stars to model PSF ({band})')
        plt.tight_layout()
        
    epsf_builder = EPSFBuilder(oversampling=oversampling, maxiters=maxiters,
                               progress_bar=False,
                               smoothing_kernel=smoothing_kernel)  
    epsf, fitted_stars = epsf_builder(stars)

    # if show:    
    #     fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize,
    #                            squeeze=True)
    #     axs = axs.ravel()
    #     for i in range(len(fitted_stars)):
    #         norm = simple_norm(fitted_stars[i], 'log', percent=99.)
    #         axs[i].imshow(fitted_stars[i], norm=norm, origin='lower', cmap='viridis')
    
    if show:
        plt.figure()
        norm = simple_norm(epsf.data, 'log', percent=99.)
        plt.imshow(epsf.data, norm=norm, origin='lower', cmap='viridis')
        plt.title(f'Effective PSF ({band})')
        plt.colorbar()
        
    return epsf

epsf = get_epsf(stars, 'g')

#%%

def get_psfphot(data, mask, psf_model, fwhm, psf_size, sigma_thres=3.,
                peakmax=None, psf_class='basic'):
    bkgrms = MADStdBackgroundRMS()
    std = bkgrms(data)
    
    threshold = sigma_thres * std
    crit_separation = 2.0 * fwhm
    
    daogroup = DAOGroup(crit_separation)
    # daogroup = DBSCANGroup(crit_separation)
    mmm_bkg = MMMBackground()
    # fitter = LMLSQFitter()
    fitter = LevMarLSQFitter()
    
    fsize = int(np.ceil(psf_size))
    fitshape = fsize if fsize%2 == 1 else fsize +1
    
    if psf_class == 'basic':
        iraffind = IRAFStarFinder(threshold=threshold,
                                  fwhm=fwhm, #sigma_psf * gaussian_sigma_to_fwhm,
                                  minsep_fwhm=0.01, roundhi=5.0, roundlo=-5.0,
                                  sharplo=0.0, sharphi=2.0, peakmax=peakmax)
        stars = iraffind(data, mask=mask)
        
        photometry = BasicPSFPhotometry(finder=iraffind,
                                        group_maker=daogroup,
                                        bkg_estimator=mmm_bkg,
                                        psf_model=psf_model,
                                        fitter=fitter,
                                        aperture_radius=1.5*fwhm,
                                        # niters=1,
                                        fitshape=fitshape)
    elif psf_class == 'daophot':
        daofind = DAOStarFinder(threshold=threshold,
                                fwhm=fwhm, #sigma_psf * gaussian_sigma_to_fwhm,
                                roundhi=5.0, roundlo=-5.0,
                                sharplo=0.0, sharphi=2.0, peakmax=peakmax)
        stars = daofind(data, mask=mask)
    
        photometry = DAOPhotPSFPhotometry(crit_separation, threshold=threshold,
                                           fwhm=fwhm, psf_model=psf_model,
                                           fitshape=fitshape,
                                           aperture_radius=1.5*fwhm,
                                           roundhi=5.0, roundlo=-5.0,
                                           sharplo=0.0, sharphi=2.0)
    
    stars.rename_columns(['xcentroid','ycentroid','flux'],
                         ['x_0', 'y_0', 'flux_0'])
    
    start = time()
    print('start time -', ctime(start))
    
    result_tab = photometry.do_photometry(image=data, mask=mask,
                                          init_guesses=stars,
                                          progress_bar=True)
    residual_image = photometry.get_residual_image()
    
    finish = time()
    print('finish time -', ctime(finish))
    print(f'elapsed time - {(finish-start)/60:.2f} min')
    
    # result_tab_dao = daophot(image=data, progress_bar=True)
    # residual_image_dao = daophot.get_residual_image()
    
    # after_dao = time()
    # print('after dao -', ctime(after_dao))
    # print(f'elapsed time - {(after_dao-after_basic)/60:.2f} min')
    # print(f'total elapsed time - {(after_dao-start)/60:.2f} min')

    shape = data.shape
    result_tab = discard_stars_outside(shape, result_tab)
    # result_tab_dao = discard_stars_outside(shape, result_tab_dao)
    # result_tab_dao = result_tab_dao[result_tab_dao['flux_fit'] > 0]

    return result_tab, photometry


def discard_stars_outside(shape, result_tab):
    ny, nx = shape
    xin = np.logical_and(result_tab['x_fit'] > 0 , result_tab['x_fit'] < nx-1)
    yin = np.logical_and(result_tab['y_fit'] > 0 , result_tab['y_fit'] < ny-1)
    isin = np.logical_and(xin, yin)
    return result_tab[isin]

import sys
sys.setrecursionlimit(10**7)

fwhm = fwhm_med
sigma_psf = fwhm * gaussian_fwhm_to_sigma
# psf_model = IntegratedGaussianPRF(sigma=sigma_psf)
psf_model = epsf
psf_size = fwhm * 4
sigma_thres = 5.
peakmax=10.0

res = get_psfphot(data, mask, psf_model, fwhm, psf_size,
                  sigma_thres=sigma_thres, peakmax=peakmax, psf_class='basic')

result_tab, photometry = res
residual_image = photometry.get_residual_image()

res = get_psfphot(data, mask, psf_model, fwhm, psf_size,
                  sigma_thres=sigma_thres, peakmax=peakmax,
                  psf_class='daophot')

result_tab_dao, photometry_dao = res
residual_image_dao = photometry_dao.get_residual_image()

hdu = fits.PrimaryHDU(residual_image)
hdu_dao = fits.PrimaryHDU(residual_image_dao)
hdu.writeto(DATADIR/f'res_{band}_{sigma_thres:.1f}.fits', overwrite=True)
hdu_dao.writeto(DATADIR/f'res_{band}_{sigma_thres:.1f}_dao.fits', overwrite=True)
result_tab.writeto(DATADIR/f'result_tab_{band}_{sigma_thres:.1f}.csv', format='csv')
result_tab_dao.writeto(DATADIR/f'result_tab_{band}_{sigma_thres:.1f}_dao.csv', format='csv')
#%%
def show_psfphot_result(data, band, result_tab, residual_image,
                        result_tab_dao, residual_image_dao):
    fig, axs = plt.subplots(2, 3, figsize=(10, 10))
    interval = ZScaleInterval()
    vmin, vmax = interval.get_limits(data)
    
    zscale_imshow(axs[0][0], data, vmin, vmax)
    zscale_imshow(axs[1][0], data, vmin, vmax)
    zscale_imshow(axs[0][1], data - residual_image, vmin, vmax)
    zscale_imshow(axs[0][2], residual_image, vmin, vmax)
    zscale_imshow(axs[1][1], data - residual_image_dao, vmin, vmax)
    zscale_imshow(axs[1][2], residual_image_dao, vmin, vmax)    
    axs[0][0].plot(result_tab['x_fit'], result_tab['y_fit'], marker='+', c='r', ms=1, ls='None')
    axs[1][0].plot(result_tab_dao['x_fit'], result_tab_dao['y_fit'], marker='+', c='r', ms=1, ls='None')
    axs[0][2].text(0.95, 0.05, f'{band}-band', ha='right', va='bottom',
                   c='w', fontweight='bold', transform=axs[0][2].transAxes)
    axs[1][2].text(0.95, 0.05, f'{band}-band', ha='right', va='bottom',
                   c='w', fontweight='bold', transform=axs[1][2].transAxes)
    axs[0][2].text(0.95, 0.01, 'Iraf star finder', ha='right', va='bottom',
                   c='w', fontweight='bold', transform=axs[0][2].transAxes)
    axs[1][2].text(0.95, 0.01, 'DAO star finder', ha='right', va='bottom',
                   c='w', fontweight='bold', transform=axs[1][2].transAxes)
    axs[0][0].set_title('Original')
    axs[0][1].set_title('Model')
    axs[0][2].set_title('Subtracted')
    axs = axs.ravel()
    for ax in axs:
        ax.axis('off')
    
    plt.tight_layout()
    
show_psfphot_result(data, band, result_tab, residual_image,
                        result_tab_dao, residual_image_dao)

#%%