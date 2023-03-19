#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 19 16:40:23 2023

Data makeup for TAO tutorial - PSF photometry

@author: Hyeonguk Bahk
"""

from astropy import units as u
from astropy.io import fits
from astropy.visualization import ZScaleInterval
from astropy.wcs import WCS
from matplotlib import pyplot as plt

from exgalcosutils.from_legacy_survey import dr2_rgb, dr10_griz_rgb

# plt.rcParams["figure.figsize"] = (7, 5)
plt.rcParams["font.family"] = 'Times New Roman'
plt.rcParams["font.size"] = 15
plt.rcParams["text.usetex"] = False
plt.rcParams["mathtext.fontset"] = 'cm'
#%%
cra, cdec = 92.8804, -69.1214

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

ax0.text(0.98, 0.02, f'$grz$ false color', transform=ax0.transAxes,
        va='bottom', ha='right', c='k', bbox=bbox)
fig.suptitle('NGC 2210 - DESI Legacy Survey $grz$ image')
