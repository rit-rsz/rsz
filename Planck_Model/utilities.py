################################################################################
# NAME : utilities
# DATE STARTED : July 9, 2020
# AUTHORS : Benjamin Vaughan
# PURPOSE : File to hold utility functions
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
from astropy.io import fits
import numpy as np
from scipy.interpolate import griddata
from math import *


def make_header(pixsize, naxis, ra, dec):
    '''
    Purpose : This function makes a fits header with all the astrometry you need to
    do transformations from ra/dec -> pixel and pixel -> ra/dec
    Inputs  : pixsize - (float) the pixel size for your map
              naxis   - (int) the size of your map (this function assumes a square map)
              ra      - (float) right ascension in degrees of the origin
              dec     - (float) declination in degrees of the origin
    Outputs : Header - a fits header object containing astrometry data
    '''
    pixsize = pixsize / 3600
    cd1_1 = -pixsize
    cd1_2 = 0
    cd2_1 = 0
    cd2_2 = pixsize
    crpix1 = naxis[0] / 2 + 1
    crpix2 = naxis[1] / 2 + 1
    crval1 = ra
    crval2 = dec
    header = fits.Header()
    header.set('NAXIS', int(2))
    header.set('NAXIS1', naxis[0])
    header.set('NAXIS2', naxis[1])
    header.set('CD1_1', cd1_1)
    header.set('CD1_2', cd1_2)
    header.set('CD2_1', cd2_1)
    header.set('CD2_2', cd2_2)
    header.set('CRPIX1', crpix1)
    header.set('CRPIX2', crpix2)
    header.set('CRVAL1', crval1)
    header.set('CRVAL2', crval2)
    header.set('BITPIX', -64)
    header.set('CTYPE1', 'RA---TAN')
    header.set('CTYPE2', 'DEC--TAN')
    return header
