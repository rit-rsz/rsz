
################################################################################
# NAME : get_spire_beam.py
# DATE STARTED : June 18, 2019
# AUTHORS : Dale Mercado
# PURPOSE :  This function computes the spire beam given a band, pixel size and
#            output map kernel.
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#      band (string)   = one of 'PSW', 'PMW', 'PLW' (def: PSW)
#      pixsize (float) = the pixel size in arcsec (def: 6/8.333/12)
#      npixx (int)     = the number of pixels required in the x axis
#                         This should be odd so the PSF is
#                         centered. (def: 5 FWHM rounded to odd)
#      npixy (int)     = the number of pixels required in the y axis.
#                         This should be odd so the PSF is
#                         centered. (def: npixx)
#      xcent (float)   = x pixel corresponding to center of the beam
#                         (can be fractional, note that pixel
#                         numbering starts at 0)  (def: npixx/2 using
#                         integral division, so if npixx is 31,
#                         this is 15).  This is in the
#                         non-oversampled beam.
#      ycent (float)   = y pixel corresponding to center of the beam
#                         (can be fractional).  (def: npixy/2, see
#                         note for xcent for further info)
# Optional inputs:
#      bolometer (string) = Optional argument specifying which bolometer
#                         to return beam for (i.e., band='PSW',
#                         bolometer='A11': psf for 'PSWA11').
#                         Not currently supported, but here for
#                         when we have bolometer specific psfs in
#                         the future.
#      fwhm (float)     = The FWHM of the beam, in arcsec.
#                          Normally this is determined by band.
#      oversamp         = Amount to oversample pixels by before
#                          convolving with pixel function.
#                          Should be an odd integer (Def: 7)
#
#
# OUTPUTS :
#           beamkern (float) = array size npixx x npixy containing
# REVISION HISTORY :
################################################################################
import scipy.io
import numpy as np
# from config import * #(this line will give me access to all directory variables)
import matplotlib.pyplot as plt
import math
from astropy.io import fits
import os
import pyfits
from astropy.modeling import models, fitting

# !!! This file is currently in an unworking state
# !!! The Gaussian Model needs to be fixed before this code can work


def get_spire_beam(band = '', pixsize = 0, npixx, npixy, xcent, ycent,\
                        bolometer = 0, FWHM=fwhm,\
                        norm = 0b, OVERSAMP=oversamp, verbose = 1, \
                        FACTOR=factor):

    if len(bolometer) != 0:
        if verbose:
            print('PSFs for specific bolometer not yet supported -- ignoring')

#   Check if we have been given a band
#   If not make an assumption
    if len(band) == 0:
        if verbose:
            print('Band parameter not supplied, assuming PSW')
        band = 'PSW'
#   Check if we have been givin a pixel size
    if len(pixsize) == 0:
        band = upper(band)
        if band == 'PSW':
            pixsize = 6
        if band == 'PMW':
            pixsize = 8 + (1/3)
        if band == 'PLW':
            pixsize = 12
        else:
            print('Unknown band'+band)
        if verbose:
            print(pixsize, FORMAT = '("pixsize parameter not '+\
                   'supplied, assuming ",F0.4," arcsec")')

    if len(fwhm) == 0:
        beamFWHM = get_spire_beam_fwhm(band)
    else:
        beamFWHM = fwhm

    if beamFWHM < 0:
        print('Invalid Beam FWHM value'+ str(beamFWHM))

#   check if we've been given the map size, if not assume something
#   npixx/npixy will be the final number of pixels

    if len(npix) == 0:
        npixx = round(beamFWHM * (5 / pixsize))
        if npixx % 2 != 1:
            npixx +=1
            if verbose:
                print('npixx not supplied, using ",I0"')
#   if no y size then assume same as x
    if len(npixy) == 0:
        npixy = npixx

# Not sure if this needs to be done in python
  # ;; make sure that these have been cast from a float properly or we get errors
    npixx = ceil(npixx)
    npixy = ceil(npixy)

    if npixx % 2 != 1 and verbose:
        print('WARNING: npixx not odd, so PSF will not be centered')
    if npixy % 2 != 1 and verbose:
        print('WARNING: npixy not odd, so psf will not be centered')

#   Now to deal with oversampling
    if len(oversamp) == 0:
        ioversamp = 7
    else:
        ioversamp = round(oversamp) #in case user provides float
        if ioversamp % 2 !=:
            print('Oversamp must be an odd interger!')

    x_gen = npixx * ioversamp
    y_gen = npixy * ioversamp
    gen_pixsize = np.float64(pixsize) / ioversamp

#   Check if we have been givin the center, if not assume middle
    if len(xcent) == 0:
        ixcent = x_gen / 2
        if verbose:
            print('xcent parameter not supplied, assuming array center')
        else:
            ixcent = xcent * ioversamp
            if ioversamp > 1:
                ixcent += ioversamp / 2

    if len(ycent) == 0:
        iycent = ygen / 2
        if verbose:
            print('ycent parameter not supplied, assuming array center')
        else:
            iycent = ycent * ioversamp
            if ioversamp > 1:
                iycent += ioversamp / 2
#   Normalize FWHM to pixels
    beamFWHM /= gen_pixsize

#   if we want this normalized than call with
    if len(factor) > 0 :
#       1D BEAM, this also now uses a scrpit called PSF_GAUSSIAN from astrolib
        beamkern = gausian


   # This is the idl args
   # psf = psf_Gaussian( NPIXEL=, FWHM= , CENTROID =
   #                  [ /DOUBLE, /NORMALIZE, ST_DEV=, NDIMEN= ] )


# Fit the data using a Gaussian
# It appears that the current Gaussian fitting paramters do not match up nicely with one that
# are used within the idl script
# g_init = models.Gaussian2D(amplitude=1, x_mean=0,\
#                            y_mean=0, x_stddev=None, y_stddev=None,\
#                            theta=None, cov_matrix=None, **kwargs)[source]¶
# fit_g = fitting.LevMarLSQFitter()
# g = fit_g(g_init, x, y)
#   Remove floating point errors???
#   If this necessary for python
    else:
        beamkerm = gaussian

        if ioversamp > 1:
            # There doesnt seem to be a python equivalent
            beamkern = rebin()


    return beamkern



if __name__ == '__main__':
    get_spire_beam('rxj1347',1)
