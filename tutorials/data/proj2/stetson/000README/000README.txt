This directory contains homegeneous photometry on Arlo Landolt's photometric
system for well-observed science targets.  The protocol follows that of the 
article

"Homogeneous photometry VII. Globular clusters in the Gaia era"

by P. B. Stetson, E. Pancino, A. Zocchi, N. Sanna, and M. Monelli, 2019 MNRAS,
485, 342

except that many of these targets are missing data in the U and/or R photometric
bandpasses.

================================================================================

The file targets.dat is an ASCII text file listing the fields that I intend to
put into this site.  

The columns are:

  mch: the total number of CCD images that I have reduced for this field

  max: the maximum number of CCD images that any single star appears in

  nmg: the total number of detections in the current catalog for the field

  fet: the number of stars that I have selected by hand as candidates for
       local photometric standards

  nU nB nV nR nI:  the actual number of stars in the field that currently meet the 
  acceptance criteria as a photometric standard, in the U, B, V, R, and I bandpasses.  The
  acceptance criteria are:
   
   1.  No fewer than five observations obtained under photometric conditions;
   2.  A standard error of the mean magnitude not greater than 0.02 mag;
   3.  Evidence of intrinsic variation in excess of readout noise, photon
       statistics, PSF fitting errors, flat-fielding errors, etc., no greater
       than 0.05 mag (root-mean-square) considering all available bandpasses.

  hU hB vV hR hI:  the relative quality of the standard stars in each field in each bandpass:
  h* = X means that in that bandpass there are at least X standard stars that have been
  observed on at least X photometric occasions.  If you want the most reliable standards,
  choose the largest values of h.

The numbers in this file are approximate, intended for guidance only.  In
particular, the number of stars in the field is almost certainly larger than the
number that can be calibrated photometrically.  If the creation dates of the
files in the folder for your particular target are more recent thant the
creation date of this file, those numbers supersede these.  If you see your
favorite field in this file and the data are not yet in this archive, please let
me know and I will deliver the data promptly.

--------------------------------------------------------------------------------

For each target there are five or six files:  

<TargetName>.pho contains calibrated photometric data for stars that I have
hand-selected to serve as local photometric standards.

<TargetName>.pos contains absolute RA and Dec and relative (x,y) positions for
the same local standards.

<TargetName>.dat contains the complete catalog of photometry and astrometry for all the
stars in the field that could be calibrated.  Columns should be self-evident.

<TargetName>.txt contains a summary of the available CCD images and credits for
the data.  The time-stamp on the file provides the vintage of the calibration.  

<TargetName>.fits is a FITS-format stacked image of each field.

<TargetName>j.fits will be present for many fields.  It is a FITS-format stacked
image of the field AFTER all fitted detections have been subtracted. 
Examination of this image will show you which stars have been missed in the
catalog.  A poorly subtracted detection may be recognized as a galaxy or image
flaw of some sort.  

Note that I am gradually populating the archive with .pos and .pho files; if you
need these and they are not present, send me an email and I will provide them.  

--------------------------------------------------------------------------------

In the catalog, the units of x and y are in arcseconds relative to an
arbitrarily chosen (0,0) for each target.  x increases east and y increases
north.  Right Ascension and Declination (2000) on the Gaia system are also
provided at the end of each line.  I have started to replace the Gaia DR2
coordinate system with Gaia EDR3; which was used is indicated at the end of
the text file provided for each field, but the difference is unimportant.

Quality indices chi and sharp are provided, as well as a Welch/Stetson
variability index and its weight for each star.

Since coordinates (0,0) in the catalog do not correspond to the corner of the
accompanying image, a transformation is required to convert from catalog
coordinates to position in the image:

     Coordinates in image pixels (X,Y) compared to catalog coordinates in
     arcseconds (x,y): 

     X = scale*(x - xo - 0.5) + 0.5
     Y = scale*(y - yo - 0.5) + 0.5

     x = (X - 0.5)/scale + xo + 0.5
     y = (Y - 0.5)/scale + yo + 0.5

     Values of xo, yo, and scale are in the FITS image header keyword OFFSCALE;
     the scale is always an integer number of pixels per arcsecond.  X increases
     east and Y increases north.  

As new images are acquired and calibrated, updated versions of these files will
replace the older ones.  If you want to watch the evolution of the photometry
for the field, keep back copies of the files you download.

================================================================================

                                 W A R N I N G !

Please be skeptical of the photometry for the very brightest stars in each
field.  Saturation effects can be subtle, and they differ from one CCD chip to
another.  They are hard to account for.  As a rule of thumb, if your science
depends critically on the absolute photometry for any individual star, and that
star is brighter than the brightest few stars in the .pho file for the field,
seek independent corroboration of the photometry wherever you can find it.

On a related note, there is always a systematic tendency for the very faintest
stars in a field to be meassured too bright.  (Luminosity functions increase
faintward, and photometric errors also increase faintward.  As a result, in the
last magnitude bin it is always easier for stars to be scattered in from below
than from above.)  If your science needs accurate photometry for the faintest
stars in my files, you would be well advised to go out and get much deeper
images of your own.  You could then use my stars (the ones that AREN'T the
faintest) to calibrate your new images.  

================================================================================

Questions, complaints, and special requests should be sent to the me:

peter.stetson@nrc-cnrc.gc.ca

Thank you.

PBS

This research used the facilities of the Canadian Astronomy Data Centre operated
by the National Research Council of Canada with the support of the Canadian
Space Agency.  

This research is based in part on data obtained from the ESO Science Archive Facility.

This research uses data obtained from the Italian Center for Astronomical Archives.

This research makes use of data obtained from the Isaac Newton Group Archive
which is maintained as part of the CASU Astronomical Data Centre at the
Institute of Astronomy, Cambridge.  

This research uses data distributed by the NOAO Science Archive. NOAO is
operated by the Association of Universities for Research in Astronomy (AURA),
Inc. under a cooperative agreement with the National Science Foundation.  

This research uses data distributed by the AAT Data Archive.

This research uses data distributed by the SMOKA Archive of NAOJ Telescopes.
