# -*- coding: utf-8 -*-
"""
@Author: Hyeonguk Bahk
@Email: bahkhyeonguk@gmail.com
@Date: 2024-03-28

@Filename: preproc.py
@Brief: This file contains the preprocessing functions for TAO tutorials.
@Details: The functions in this file are responsible for preprocessing data
before further analysis.
"""

import numpy as np
from pathlib import Path
import re
import time
import warnings

from astropy.io import fits
from astropy.table import Table
from astropy.visualization import ZScaleInterval
from astropy.nddata import CCDData
from astropy.stats import mad_std
import astroalign as aa
from ccdproc import combine
from matplotlib import pyplot as plt


def make_summary_table_sao_port1(rawdir, suffix='.fit'):
    """Make a summary table of the raw data in the directory. This function is
    specifically designed for the data taken with the SAO 1m Telscope Port1.

    Args:
        rawdir (pathlib.Path): The directory containing the raw data.
        suffix (str, optional): The suffix of the raw data files. Defaults to '.fit'.

    Returns:
        stab (astropy.table.Table): The summary table of the raw data in the directory.
        
    """
    # making a summary table
    summary = []
    for f in rawdir.glob('*'+suffix):
        hdr = fits.getheader(f)
        
        # getting the filter information from the header
        filt = 'unknown'
        if 'filter' in hdr:
            filt = hdr['filter']
        elif 'FILTER' in hdr:
            filt = hdr['FILTER']
        else:
            # checking for the filter in the filename with regex, if not in header
            mtch = re.search(r'[UBVRIgriz].fit$', f.name)
            mtch_Ha = re.search(r'Ha.fit$', f.name)
            if mtch:
                filt = mtch.group(0)[0]
            elif mtch_Ha:
                filt = 'Ha'
            
        # getting the airmass X (Pickering 2002)
        try:
            alt = float(hdr['OBJCTALT'])
            X = 1/(np.sin(np.radians( alt + 244/(165 + 47*alt**1.1 ))))
        except TypeError:
            X = -1
        except Exception as e:
            print(f'Error in {f.name}: {e}')
            X = -1
        
        # appending the data to the summary list
        summary.append([f.name, hdr['DATE-OBS'], hdr['OBJECT'],
                        hdr['OBJCTRA'], hdr['OBJCTDEC'],
                        hdr['IMAGETYP'], hdr['EXPTIME'], X, filt])
        
    # creating the summary table
    stab = Table(rows=summary,
                 names=['filename', 'date-obs', 'object', 'ra', 'dec', 'imagetyp',
                        'exptime', 'airmass', 'filter'],
                 dtype=['U50', 'U50', 'U50', 'U50', 'U50', 'U50', 'f8', 'f8', 'U50'])
    return stab


def combine_bias(biaslist, outdir, outname, unit='adu'):
    """Combine the bias frames in the directory to make a master bias frame.

    Args:
        biaslist (list[Path]): The list of paths to the bias frames.
        outdir (Path): The directory to save the combined bias frame.
        outname (str): The name of the combined bias frame.
        unit (str, optional): The unit of the CCD data to pass to the CCDData object.
            Defaults to 'adu'.

    Returns:
        CCDData: The combined bias frame.

    """
    # reading the bias frames
    bias = [CCDData.read(f, unit=unit) for f in biaslist]

    # combining the bias frames
    mbias = combine(bias, method='average', sigma_clip=True,
                    sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                    sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std)

    # calculating the readout noise
    if 'EGAIN' in bias[0].header:
        gain = bias[0].header['EGAIN'] # gain in e-/ADU
    elif 'GAIN' in bias[0].header:
        gain = bias[0].header['GAIN'] # gain in e-/ADU
    else:
        raise ValueError('Gain not found in the header.')
    bias1 = bias[0].data.astype('int16')
    bias2 = bias[1].data.astype('int16')
    rdnoise = np.std(bias1 - bias2) / np.sqrt(2) * gain  # in e-

    # writing the readout noise to the header
    mbias.header['RDNOISE'] = (rdnoise, 'Readout noise in e-')

    # saving the combined bias frame
    try:
        mbias.write(outdir / outname, overwrite=True)
    except Exception as e:
        print(f'Error saving the combined bias frame: {e}')

    return mbias


def make_master_dark(dark_list, outdir, outname=None, mbias=None, verbose=True):
    """Combine and process a list of dark frames to create a master dark frame.

    Args:
        dark_list (list[Path]): The list of paths to the dark frames.
        outdir (Path): The directory to save the master dark frame.
        outname (str): The name of the master dark frame. If None, the default name
            will be `MDark{exptime}.fits`. Defaults to None.
        mbias (Path, optional): The master bias frame to subtract from the dark frames.
            Defaults to None.
        verbose (bool, optional): Whether to print verbose output. Defaults to True.

    Returns:
        CCDData: The master dark frame.

    """
    # - checking the basic info
    # - Please check the consistency in observation dates and exposure times.
    if verbose:
        for i, dark_file in enumerate(dark_list):
            dark_hdr = fits.getheader(dark_file)
            print(f"\nDark frame {i+1:d}")
            for key in ['DATE-OBS', 'EXPTIME']:
                print(f"  {key} = {dark_hdr[key]}")
    
    if mbias is None:
        mbias_data = 0.
    else:
        mbias_data = fits.getdata(mbias)
        
    # - stacking dark frames
    dark_stack = []
    start = time.time()
    for dark_file in dark_list:
        dark_data, dark_hdr = fits.getdata(dark_file, header=True)
        dark_bn = (dark_data - mbias_data)
        dark = CCDData(data=dark_bn, header=dark_hdr, unit='adu')    
        dark_stack.append(dark)
    finish = time.time()
    print(f"\nReading {len(dark_list)} dark frames took {finish-start:.2f} sec")
        
    # - combine with sigma clipping
    start = time.time()
    mdark = combine(dark_stack, sigma_clip=True,
                    sigma_clip_high_thresh=3, sigma_clip_low_thresh=3)
    finish = time.time()
    print(f"Combining {len(dark_list)} dark frames took {finish-start:.2f} sec")
    
    # - correcting the negative values
    mdark.data[mdark.data < 0.] = 0.
    
    # - save the master dark as fits file
    start = time.time()
    dark_hdr['NFRAMES'] = len(dark_list)  # recording # of dark frames combined
    if outname is None:
        outname = f"MDark{dark_hdr['EXPTIME']:.0f}.fits"
    fits.writeto(outdir/outname, mdark, dark_hdr, overwrite=True)
    finish = time.time()
    print(f"Writing the master dark took {finish-start:.2f} sec")
        
    return mdark


def make_master_flat(flat_list, outdir, outname=None, mbias=None, mdark=None,
                     unit='adu', verbose=True, filter_key='FILTER'):
    """Combine and process a list of flat frames to create a master flat frame.

    Args:
        flat_list (list[Path]): The list of paths to the flat frames.
        outdir (Path): The directory to save the master flat frame.
        outname (str): The name of the master flat frame. If None, the default name
            will be `MFlat{filter}.fits`. Defaults to None.
        mbias (Path, optional): The master bias frame to subtract from the flat frames.
            Defaults to None.
        mdark (Path, optional): The master dark frame to subtract from the flat frames.
            Defaults to None.
        unit (str, optional): The unit of the CCD data to pass to the CCDData object.
            Defaults to 'adu'.
        verbose (bool, optional): Whether to print verbose output. Defaults to True.
        filter_key (str, optional): The header keyword for the filter information.
            This is used when the verbose option is True or when the output name is
            not provided. Defaults to 'FILTER'.

    Returns:
        CCDData: The master flat frame.

    """
    if verbose:
        # checking the basic info: check dates, exposure times, and filters
        for i, flat_file in enumerate(flat_list):
            flat_hdr = fits.getheader(flat_file)
            print(f"\nFlat frame {i+1:d}")
            for key in ['DATE-OBS', 'EXPTIME', filter_key]:
                print(f"  {key} = {flat_hdr[key]}")
    
    # reading the master bias and dark frames
    mbias_data = 0. if mbias is None else fits.getdata(mbias)
    mdark_data = 0. if mdark is None else fits.getdata(mdark)
    
    # stacking flat frames
    flat_stack = []   
    for i, flat_file in enumerate(flat_list):
        flat_data, flat_hdr = fits.getdata(flat_file, header=True)  
        # bias and dark subtraction
        flat_bd = (flat_data - mbias_data - mdark_data)
        # flat scaling (with relative sensitivity=1 at the maximum)
        flat_bdn = flat_bd / np.median(flat_bd)
        flat_stack.append(CCDData(data=flat_bdn, unit=unit))
            
    # sigma clipping
    if verbose:
        print(f"\nCombining flat frames with sigma clipping... {len(flat_list)} frames")
    mflat = combine(flat_stack, sigma_clip=True,
                    sigma_clip_low_thresh=3, sigma_clip_high_thresh=3,
                    sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std)
    
    # save the master flat as fits file
    # Caveat: this assumes that the flat frames are taken with the same
    #         filter, and they have the same exposure time. If not, you may
    #         need to modify the filename.
    flat_hdr['NFRAMES'] = len(flat_list)
    if outname is None:
        filter_now = flat_hdr[filter_key]     # specifying current filter
        outname = f"MFlat{filter_now}.fits"
    fits.writeto(outdir / outname, mflat.data, header=flat_hdr, overwrite=True)
        
    return mflat

def preproc(sci_list, outdir, outname=None, mbias=None, mdark=None, mflat=None,
            rdnoise=None, verbose=True):
    """ Preprocesses a list of science images by performing bias subtraction, dark
    subtraction, and flat fielding.

    Args:
        sci_list (list[Path]): List of paths to science images.
        outdir (Path): Output directory to save the preprocessed images.
        outname (str, optional): Output file name for the preprocessed images. If not
            provided, a default name will be prepended with 'p'. Defaults to None.
        mbias (Path, optional): Path to the master bias frame. If not provided, bias
            subtraction will not be performed.
        mdark (Path, optional): Path to the master dark frame. If not provided, dark
            subtraction will not be performed.
        mflat (Path, optional): Path to the master flat frame. If not provided, flat
            fielding will not be performed.
        rdnoise (float, optional): Readout noise in electrons. If not provided, it will
            be extracted from the master bias frame header.
        verbose (bool, optional): Whether to print progress messages. Defaults to True.

    Returns:
        None
    """
    for i in range(len(sci_list)):
        # bias subtraction, dark subtraction, and flat fielding
        sci_path = sci_list[i]
        sci_data, sci_hdr  = fits.getdata(sci_path, header=True)
        # 'int' type may cause error when calculating
        sci_data0 = sci_data.astype('float')
        
        mbias_data = 0. if mbias is None else fits.getdata(mbias)
        mdark_data = 0. if mdark is None else fits.getdata(mdark)
        mflat_data = 1. if mflat is None else fits.getdata(mflat)
        
        sci_data = sci_data0 - mbias_data    # Bias subtraction
        sci_data -= mdark_data   # Dark subtraction
        sci_data /= mflat_data    # Flat fielding
        
        # recording preprocessing history
        now = time.strftime("%Y-%m-%d %H:%M:%S (GMT%z)")
        if rdnoise is None:
            mbias_hdr = fits.getheader(mbias)
            rdnoise = mbias_hdr['RDNOISE']
        sci_hdr['RDNOISE'] = (rdnoise, 'Readout noise in e-')
        sci_hdr['HISTORY'] = 'Preprocessed at ' + now
            
        # saving preprocessed image to a fits file
        if outname is None:
            outname = f"p{sci_path.stem}.fits"
        fits.writeto(outdir / outname, sci_data, sci_hdr, overwrite=True)
        if verbose:
            print(f"Done: {sci_path.stem}")
