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

def save_header(ref_head, avg_b, avg_T, avg_tau, id, unit, cent):
    header = fits.Header()
    header.set('UNITS', unit)
    header.set('CRVAL1', cent[0])
    header.set('CRVAL2', cent[1])
    header.set('NAXIS', ref_head['NAXIS'])
    header.set('NAXIS1', ref_head['NAXIS1'])
    header.set('NAXIS2', ref_head['NAXIS2'])
    header.set('CD1_1', ref_head['CD1_1'], 'wcs pixel scales')
    header.set('CD1_2', ref_head['CD1_2'])
    header.set('CD2_1', ref_head['CD2_1'])
    header.set('CD2_2', ref_head['Cd2_2'])
    header.set('CRPIX1', ref_head['NAXIS1'] / 2, 'center pixel')
    header.set('CRPIX2', ref_head['NAXIS2'] / 2)
    header.set('CRVAL1', ref_head['CRVAL1'], 'center pixel in degrees')
    header.set('CRVAL2', ref_head['CRVAL2'])
    header.set('BITPIX', -64)
    header.set('CTYPE1', 'RA---TAN')
    header.set('CTYPE2', 'DEC--TAN')
    header.set('RADESYS', 'FK5', 'reference frame')
    header.set('AVGTAU', avg_tau, 'average value for optical depth')
    header.set('AVGTEMP', avg_T, 'average temperature')
    header.set('AVGBETA', avg_b, 'average beta')
    return header
