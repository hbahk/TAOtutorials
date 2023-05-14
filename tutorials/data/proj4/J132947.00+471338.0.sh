#!/bin/bash
#
# This script will create a mosaic of SDSS images.
#
# To invoke this script do:
#
# $ /bin/bash J132947.00+471338.0.sh
#
hasPhotoObj=''
if [[ -n "${BOSS_PHOTOOBJ}" ]]; then
    if [[ -d ${BOSS_PHOTOOBJ}/frames/301 ]]; then
        echo "BOSS_PHOTOOBJ detected.  Will attempt to use files on disk."
        hasPhotoObj=True
    fi
fi
for cmd in swarp bzip2 wget; do
    hasCmd=$(which ${cmd} 2>/dev/null)
    if [[ -z "${hasCmd}" ]]; then
        echo "This script requires ${cmd}, which is not in your \$PATH."
        exit 1
    fi
done
if [[ -z "${hasPhotoObj}" ]]; then
    if [[ ! -f frame-r-003650-1-0090.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3650/1/frame-r-003650-1-0090.fits.bz2
        bzip2 -d -c -k frame-r-003650-1-0090.fits.bz2 > frame-r-003650-1-0090.fits
    fi
    if [[ ! -f frame-r-003699-5-0099.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3699/5/frame-r-003699-5-0099.fits.bz2
        bzip2 -d -c -k frame-r-003699-5-0099.fits.bz2 > frame-r-003699-5-0099.fits
    fi
    if [[ ! -f frame-r-003699-5-0100.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3699/5/frame-r-003699-5-0100.fits.bz2
        bzip2 -d -c -k frame-r-003699-5-0100.fits.bz2 > frame-r-003699-5-0100.fits
    fi
    if [[ ! -f frame-r-003699-5-0101.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3699/5/frame-r-003699-5-0101.fits.bz2
        bzip2 -d -c -k frame-r-003699-5-0101.fits.bz2 > frame-r-003699-5-0101.fits
    fi
    if [[ ! -f frame-r-003716-5-0117.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3716/5/frame-r-003716-5-0117.fits.bz2
        bzip2 -d -c -k frame-r-003716-5-0117.fits.bz2 > frame-r-003716-5-0117.fits
    fi
    if [[ ! -f frame-r-003716-5-0119.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3716/5/frame-r-003716-5-0119.fits.bz2
        bzip2 -d -c -k frame-r-003716-5-0119.fits.bz2 > frame-r-003716-5-0119.fits
    fi
    if [[ ! -f frame-r-003716-6-0117.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3716/6/frame-r-003716-6-0117.fits.bz2
        bzip2 -d -c -k frame-r-003716-6-0117.fits.bz2 > frame-r-003716-6-0117.fits
    fi
    if [[ ! -f frame-r-003716-6-0119.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3716/6/frame-r-003716-6-0119.fits.bz2
        bzip2 -d -c -k frame-r-003716-6-0119.fits.bz2 > frame-r-003716-6-0119.fits
    fi
    if [[ ! -f frame-r-002964-1-0402.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/2964/1/frame-r-002964-1-0402.fits.bz2
        bzip2 -d -c -k frame-r-002964-1-0402.fits.bz2 > frame-r-002964-1-0402.fits
    fi
    if [[ ! -f frame-r-003699-6-0098.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3699/6/frame-r-003699-6-0098.fits.bz2
        bzip2 -d -c -k frame-r-003699-6-0098.fits.bz2 > frame-r-003699-6-0098.fits
    fi
    if [[ ! -f frame-r-003716-5-0114.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3716/5/frame-r-003716-5-0114.fits.bz2
        bzip2 -d -c -k frame-r-003716-5-0114.fits.bz2 > frame-r-003716-5-0114.fits
    fi
    if [[ ! -f frame-r-003716-6-0114.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3716/6/frame-r-003716-6-0114.fits.bz2
        bzip2 -d -c -k frame-r-003716-6-0114.fits.bz2 > frame-r-003716-6-0114.fits
    fi
    if [[ ! -f frame-r-003716-6-0115.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3716/6/frame-r-003716-6-0115.fits.bz2
        bzip2 -d -c -k frame-r-003716-6-0115.fits.bz2 > frame-r-003716-6-0115.fits
    fi
    if [[ ! -f frame-r-003650-1-0086.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3650/1/frame-r-003650-1-0086.fits.bz2
        bzip2 -d -c -k frame-r-003650-1-0086.fits.bz2 > frame-r-003650-1-0086.fits
    fi
    if [[ ! -f frame-r-003699-6-0103.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3699/6/frame-r-003699-6-0103.fits.bz2
        bzip2 -d -c -k frame-r-003699-6-0103.fits.bz2 > frame-r-003699-6-0103.fits
    fi
    if [[ ! -f frame-r-003716-5-0120.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3716/5/frame-r-003716-5-0120.fits.bz2
        bzip2 -d -c -k frame-r-003716-5-0120.fits.bz2 > frame-r-003716-5-0120.fits
    fi
    if [[ ! -f frame-r-003716-6-0120.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3716/6/frame-r-003716-6-0120.fits.bz2
        bzip2 -d -c -k frame-r-003716-6-0120.fits.bz2 > frame-r-003716-6-0120.fits
    fi
    if [[ ! -f frame-r-003699-6-0099.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3699/6/frame-r-003699-6-0099.fits.bz2
        bzip2 -d -c -k frame-r-003699-6-0099.fits.bz2 > frame-r-003699-6-0099.fits
    fi
    if [[ ! -f frame-r-003699-6-0100.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3699/6/frame-r-003699-6-0100.fits.bz2
        bzip2 -d -c -k frame-r-003699-6-0100.fits.bz2 > frame-r-003699-6-0100.fits
    fi
    if [[ ! -f frame-r-003699-6-0101.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3699/6/frame-r-003699-6-0101.fits.bz2
        bzip2 -d -c -k frame-r-003699-6-0101.fits.bz2 > frame-r-003699-6-0101.fits
    fi
    if [[ ! -f frame-r-003699-6-0102.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3699/6/frame-r-003699-6-0102.fits.bz2
        bzip2 -d -c -k frame-r-003699-6-0102.fits.bz2 > frame-r-003699-6-0102.fits
    fi
    if [[ ! -f frame-r-003716-5-0116.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3716/5/frame-r-003716-5-0116.fits.bz2
        bzip2 -d -c -k frame-r-003716-5-0116.fits.bz2 > frame-r-003716-5-0116.fits
    fi
    if [[ ! -f frame-r-003716-5-0118.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3716/5/frame-r-003716-5-0118.fits.bz2
        bzip2 -d -c -k frame-r-003716-5-0118.fits.bz2 > frame-r-003716-5-0118.fits
    fi
    if [[ ! -f frame-r-003716-6-0116.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3716/6/frame-r-003716-6-0116.fits.bz2
        bzip2 -d -c -k frame-r-003716-6-0116.fits.bz2 > frame-r-003716-6-0116.fits
    fi
    if [[ ! -f frame-r-003650-1-0087.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3650/1/frame-r-003650-1-0087.fits.bz2
        bzip2 -d -c -k frame-r-003650-1-0087.fits.bz2 > frame-r-003650-1-0087.fits
    fi
    if [[ ! -f frame-r-003650-1-0088.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3650/1/frame-r-003650-1-0088.fits.bz2
        bzip2 -d -c -k frame-r-003650-1-0088.fits.bz2 > frame-r-003650-1-0088.fits
    fi
    if [[ ! -f frame-r-003716-6-0118.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3716/6/frame-r-003716-6-0118.fits.bz2
        bzip2 -d -c -k frame-r-003716-6-0118.fits.bz2 > frame-r-003716-6-0118.fits
    fi
    if [[ ! -f frame-r-002964-1-0400.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/2964/1/frame-r-002964-1-0400.fits.bz2
        bzip2 -d -c -k frame-r-002964-1-0400.fits.bz2 > frame-r-002964-1-0400.fits
    fi
    if [[ ! -f frame-r-002964-1-0401.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/2964/1/frame-r-002964-1-0401.fits.bz2
        bzip2 -d -c -k frame-r-002964-1-0401.fits.bz2 > frame-r-002964-1-0401.fits
    fi
    if [[ ! -f frame-r-003650-1-0089.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3650/1/frame-r-003650-1-0089.fits.bz2
        bzip2 -d -c -k frame-r-003650-1-0089.fits.bz2 > frame-r-003650-1-0089.fits
    fi
    if [[ ! -f frame-r-003699-6-0097.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3699/6/frame-r-003699-6-0097.fits.bz2
        bzip2 -d -c -k frame-r-003699-6-0097.fits.bz2 > frame-r-003699-6-0097.fits
    fi
    if [[ ! -f frame-r-003716-5-0115.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3716/5/frame-r-003716-5-0115.fits.bz2
        bzip2 -d -c -k frame-r-003716-5-0115.fits.bz2 > frame-r-003716-5-0115.fits
    fi
    if [[ ! -f frame-r-002964-1-0399.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/2964/1/frame-r-002964-1-0399.fits.bz2
        bzip2 -d -c -k frame-r-002964-1-0399.fits.bz2 > frame-r-002964-1-0399.fits
    fi
    if [[ ! -f frame-r-003650-1-0085.fits ]]; then
        wget -q https://data.sdss.org/sas/dr12/env/BOSS_PHOTOOBJ/frames/301/3650/1/frame-r-003650-1-0085.fits.bz2
        bzip2 -d -c -k frame-r-003650-1-0085.fits.bz2 > frame-r-003650-1-0085.fits
    fi
else
    if [[ ! -f frame-r-003650-1-0090.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3650/1/frame-r-003650-1-0090.fits.bz2 > frame-r-003650-1-0090.fits
    fi
    if [[ ! -f frame-r-003699-5-0099.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3699/5/frame-r-003699-5-0099.fits.bz2 > frame-r-003699-5-0099.fits
    fi
    if [[ ! -f frame-r-003699-5-0100.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3699/5/frame-r-003699-5-0100.fits.bz2 > frame-r-003699-5-0100.fits
    fi
    if [[ ! -f frame-r-003699-5-0101.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3699/5/frame-r-003699-5-0101.fits.bz2 > frame-r-003699-5-0101.fits
    fi
    if [[ ! -f frame-r-003716-5-0117.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3716/5/frame-r-003716-5-0117.fits.bz2 > frame-r-003716-5-0117.fits
    fi
    if [[ ! -f frame-r-003716-5-0119.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3716/5/frame-r-003716-5-0119.fits.bz2 > frame-r-003716-5-0119.fits
    fi
    if [[ ! -f frame-r-003716-6-0117.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3716/6/frame-r-003716-6-0117.fits.bz2 > frame-r-003716-6-0117.fits
    fi
    if [[ ! -f frame-r-003716-6-0119.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3716/6/frame-r-003716-6-0119.fits.bz2 > frame-r-003716-6-0119.fits
    fi
    if [[ ! -f frame-r-002964-1-0402.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/2964/1/frame-r-002964-1-0402.fits.bz2 > frame-r-002964-1-0402.fits
    fi
    if [[ ! -f frame-r-003699-6-0098.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3699/6/frame-r-003699-6-0098.fits.bz2 > frame-r-003699-6-0098.fits
    fi
    if [[ ! -f frame-r-003716-5-0114.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3716/5/frame-r-003716-5-0114.fits.bz2 > frame-r-003716-5-0114.fits
    fi
    if [[ ! -f frame-r-003716-6-0114.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3716/6/frame-r-003716-6-0114.fits.bz2 > frame-r-003716-6-0114.fits
    fi
    if [[ ! -f frame-r-003716-6-0115.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3716/6/frame-r-003716-6-0115.fits.bz2 > frame-r-003716-6-0115.fits
    fi
    if [[ ! -f frame-r-003650-1-0086.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3650/1/frame-r-003650-1-0086.fits.bz2 > frame-r-003650-1-0086.fits
    fi
    if [[ ! -f frame-r-003699-6-0103.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3699/6/frame-r-003699-6-0103.fits.bz2 > frame-r-003699-6-0103.fits
    fi
    if [[ ! -f frame-r-003716-5-0120.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3716/5/frame-r-003716-5-0120.fits.bz2 > frame-r-003716-5-0120.fits
    fi
    if [[ ! -f frame-r-003716-6-0120.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3716/6/frame-r-003716-6-0120.fits.bz2 > frame-r-003716-6-0120.fits
    fi
    if [[ ! -f frame-r-003699-6-0099.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3699/6/frame-r-003699-6-0099.fits.bz2 > frame-r-003699-6-0099.fits
    fi
    if [[ ! -f frame-r-003699-6-0100.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3699/6/frame-r-003699-6-0100.fits.bz2 > frame-r-003699-6-0100.fits
    fi
    if [[ ! -f frame-r-003699-6-0101.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3699/6/frame-r-003699-6-0101.fits.bz2 > frame-r-003699-6-0101.fits
    fi
    if [[ ! -f frame-r-003699-6-0102.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3699/6/frame-r-003699-6-0102.fits.bz2 > frame-r-003699-6-0102.fits
    fi
    if [[ ! -f frame-r-003716-5-0116.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3716/5/frame-r-003716-5-0116.fits.bz2 > frame-r-003716-5-0116.fits
    fi
    if [[ ! -f frame-r-003716-5-0118.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3716/5/frame-r-003716-5-0118.fits.bz2 > frame-r-003716-5-0118.fits
    fi
    if [[ ! -f frame-r-003716-6-0116.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3716/6/frame-r-003716-6-0116.fits.bz2 > frame-r-003716-6-0116.fits
    fi
    if [[ ! -f frame-r-003650-1-0087.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3650/1/frame-r-003650-1-0087.fits.bz2 > frame-r-003650-1-0087.fits
    fi
    if [[ ! -f frame-r-003650-1-0088.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3650/1/frame-r-003650-1-0088.fits.bz2 > frame-r-003650-1-0088.fits
    fi
    if [[ ! -f frame-r-003716-6-0118.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3716/6/frame-r-003716-6-0118.fits.bz2 > frame-r-003716-6-0118.fits
    fi
    if [[ ! -f frame-r-002964-1-0400.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/2964/1/frame-r-002964-1-0400.fits.bz2 > frame-r-002964-1-0400.fits
    fi
    if [[ ! -f frame-r-002964-1-0401.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/2964/1/frame-r-002964-1-0401.fits.bz2 > frame-r-002964-1-0401.fits
    fi
    if [[ ! -f frame-r-003650-1-0089.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3650/1/frame-r-003650-1-0089.fits.bz2 > frame-r-003650-1-0089.fits
    fi
    if [[ ! -f frame-r-003699-6-0097.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3699/6/frame-r-003699-6-0097.fits.bz2 > frame-r-003699-6-0097.fits
    fi
    if [[ ! -f frame-r-003716-5-0115.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3716/5/frame-r-003716-5-0115.fits.bz2 > frame-r-003716-5-0115.fits
    fi
    if [[ ! -f frame-r-002964-1-0399.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/2964/1/frame-r-002964-1-0399.fits.bz2 > frame-r-002964-1-0399.fits
    fi
    if [[ ! -f frame-r-003650-1-0085.fits ]]; then
        bzip2 -d -c -k ${BOSS_PHOTOOBJ}/frames/301/3650/1/frame-r-003650-1-0085.fits.bz2 > frame-r-003650-1-0085.fits
    fi
fi
#
# Filter = r
#
/bin/rm -f default.swarp
cat > default.swarp <<EOT
IMAGEOUT_NAME          J132947.00+471338.0-r.fits      # Output filename
WEIGHTOUT_NAME       J132947.00+471338.0-r.weight.fits # Output weight-map filename

HEADER_ONLY            N               # Only a header as an output file (Y/N)?
HEADER_SUFFIX          .head           # Filename extension for additional headers

#------------------------------- Input Weights --------------------------------

WEIGHT_TYPE            NONE            # BACKGROUND,MAP_RMS,MAP_VARIANCE
                                       # or MAP_WEIGHT
WEIGHT_SUFFIX          weight.fits     # Suffix to use for weight-maps
WEIGHT_IMAGE                           # Weightmap filename if suffix not used
                                       # (all or for each weight-map)

#------------------------------- Co-addition ----------------------------------

COMBINE                Y               # Combine resampled images (Y/N)?
COMBINE_TYPE           MEDIAN          # MEDIAN,AVERAGE,MIN,MAX,WEIGHTED,CHI2
                                       # or SUM

#-------------------------------- Astrometry ----------------------------------

CELESTIAL_TYPE         NATIVE          # NATIVE, PIXEL, EQUATORIAL,
                                       # GALACTIC,ECLIPTIC, or SUPERGALACTIC
PROJECTION_TYPE        TAN             # Any WCS projection code or NONE
PROJECTION_ERR         0.001           # Maximum projection error (in output
                                       # pixels), or 0 for no approximation
CENTER_TYPE            MANUAL          # MANUAL, ALL or MOST
CENTER                       202.4465216670,       47.2273441667 # Image Center
PIXELSCALE_TYPE        MANUAL          # MANUAL,FIT,MIN,MAX or MEDIAN
PIXEL_SCALE            0.396127  # Pixel scale
IMAGE_SIZE             2726,2726 # scale = 0.396127 arcsec/pixel

#-------------------------------- Resampling ----------------------------------

RESAMPLE               Y               # Resample input images (Y/N)?
RESAMPLE_DIR           .               # Directory path for resampled images
RESAMPLE_SUFFIX        .resamp.fits    # filename extension for resampled images

RESAMPLING_TYPE        LANCZOS3        # NEAREST,BILINEAR,LANCZOS2,LANCZOS3
                                       # or LANCZOS4 (1 per axis)
OVERSAMPLING           0               # Oversampling in each dimension
                                       # (0 = automatic)
INTERPOLATE            N               # Interpolate bad input pixels (Y/N)?
                                       # (all or for each image)

FSCALASTRO_TYPE        FIXED           # NONE,FIXED, or VARIABLE
FSCALE_KEYWORD         FLXSCALE        # FITS keyword for the multiplicative
                                       # factor applied to each input image
FSCALE_DEFAULT         1.0             # Default FSCALE value if not in header

GAIN_KEYWORD           GAIN            # FITS keyword for effect. gain (e-/ADU)
GAIN_DEFAULT           0.0             # Default gain if no FITS keyword found

#--------------------------- Background subtraction ---------------------------

SUBTRACT_BACK          N               # Subtraction sky background (Y/N)?
                                       # (all or for each image)

BACK_TYPE              AUTO            # AUTO or MANUAL
                                       # (all or for each image)
BACK_DEFAULT           0.0             # Default background value in MANUAL
                                       # (all or for each image)
BACK_SIZE              128             # Background mesh size (pixels)
                                       # (all or for each image)
BACK_FILTERSIZE        3               # Background map filter range (meshes)
                                       # (all or for each image)

#------------------------------ Memory management -----------------------------

VMEM_DIR               .               # Directory path for swap files
VMEM_MAX               2047            # Maximum amount of virtual memory (MB)
MEM_MAX                2048            # Maximum amount of usable RAM (MB)
COMBINE_BUFSIZE        1024            # Buffer size for combine (MB)

#------------------------------ Miscellaneous ---------------------------------

DELETE_TMPFILES        Y               # Delete temporary resampled FITS files
                                       # (Y/N)?
COPY_KEYWORDS          OBJECT          # List of FITS keywords to propagate
                                       # from the input to the output headers
WRITE_FILEINFO         Y               # Write information about each input
                                       # file in the output image header?
WRITE_XML              N               # Write XML file (Y/N)?
XML_NAME               swarp.xml       # Filename for XML output
VERBOSE_TYPE           QUIET           # QUIET,NORMAL or FULL

NTHREADS               0               # Number of simultaneous threads for
                                       # the SMP version of SWarp
                                       # 0 = automatic

EOT
#
swarp frame-r-003650-1-0090.fits frame-r-003699-5-0099.fits frame-r-003699-5-0100.fits frame-r-003699-5-0101.fits frame-r-003716-5-0117.fits frame-r-003716-5-0119.fits frame-r-003716-6-0117.fits frame-r-003716-6-0119.fits frame-r-002964-1-0402.fits frame-r-003699-6-0098.fits frame-r-003716-5-0114.fits frame-r-003716-6-0114.fits frame-r-003716-6-0115.fits frame-r-003650-1-0086.fits frame-r-003699-6-0103.fits frame-r-003716-5-0120.fits frame-r-003716-6-0120.fits frame-r-003699-6-0099.fits frame-r-003699-6-0100.fits frame-r-003699-6-0101.fits frame-r-003699-6-0102.fits frame-r-003716-5-0116.fits frame-r-003716-5-0118.fits frame-r-003716-6-0116.fits frame-r-003650-1-0087.fits frame-r-003650-1-0088.fits frame-r-003716-6-0118.fits frame-r-002964-1-0400.fits frame-r-002964-1-0401.fits frame-r-003650-1-0089.fits frame-r-003699-6-0097.fits frame-r-003716-5-0115.fits frame-r-002964-1-0399.fits frame-r-003650-1-0085.fits
