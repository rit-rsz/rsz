
################################################################################
# NAME : get_spire_beam.py
# DATE STARTED : June 18, 2019
# AUTHORS : Dale Mercado & Benjamin Vaugahn
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
import scipy.signal
import os
# from astropy.modeling import models, fitting
from astropy.convolution import Gaussian2DKernel
from get_spire_beam_fwhm import *


# !!! This file is currently in an unworking state
# !!! The Gaussian Model needs to be fixed before this code can work


def get_spire_beam(band = '',pixsize = 0,npixx=0, npixy=0,\
                   xcent=0, ycent=0,bolometer=0, fwhm='',\
                   norm=0, oversamp=0, verbose=1,\
                   factor=0):
    if bolometer != 0:
        if verbose:
            print('PSFs for specific bolometer not yet supported -- ignoring')

#   Check if we have been given a band
#   If not make an assumption
    if len(band) == 0:
        if verbose:
            print('Band parameter not supplied, assuming PSW')
        band = 'PSW'
#   Check if we have been givin a pixel size
    if pixsize == 0:
        # band = upper(band)
        if band == 'PSW':
            pixsize = 6
        if band == 'PMW':
            pixsize = 8 + (1/3)
        if band == 'PLW':
            pixsize = 12
        else:
            print('Unknown band'+band)
        if verbose:
            print('pixsize paramter not supplied assuming %s arcsec' % (pixsize))

    # Figure out which color this is

    if len(fwhm) == 0:
        beamFWHM = get_spire_beam_fwhm(band)
    else:
        beamFWHM = fwhm

    if beamFWHM < 0:
        print('Invalid Beam FWHM value'+ str(beamFWHM))

#   check if we've been given the map size, if not assume something
#   npixx/npixy will be the final number of pixels

    if npixx == 0:
        npixx = round(beamFWHM * 5 / (pixsize))
        if npixx % 2 != 1:
            npixx +=1
            if verbose:
                print('npixx not supplied, using ",I0"')
    #If no y size then assume same as x
    if npixy == 0:
        npixy = npixx

    #Not sure if this needs to be done in python
    #Make sure that these have been cast from a float properly or we get errors
    npixx = math.ceil(npixx)
    npixy = math.ceil(npixy)

    if npixx % 2 != 1 and verbose:
        print('WARNING: npixx not odd, so PSF will not be centered')
    if npixy % 2 != 1 and verbose:
        print('WARNING: npixy not odd, so psf will not be centered')

#   Now to deal with oversampling
    if oversamp == 0:
        ioversamp = 7
    else:
        ioversamp = round(oversamp) #in case user provides float
        if ioversamp % 2 != 1:
            print('Oversamp must be an odd interger!')

    x_gen = npixx * ioversamp
    y_gen = npixy * ioversamp
    gen_pixsize = np.float64(pixsize) / ioversamp

#   Check if we have been givin the center, if not assume middle
    if xcent == 0:
        ixcent = x_gen / 2
        if verbose:
            print('xcent parameter not supplied, assuming array center')
    else:
        ixcent = xcent * ioversamp
        if ioversamp > 1:
            ixcent = ixcent + ioversamp / 2

    if ycent == 0:
        iycent = y_gen / 2
        if verbose:
            print('ycent parameter not supplied, assuming array center')
    else:
        iycent = ycent * ioversamp
        if ioversamp > 1:
            iycent = iycent +ioversamp / 2
#   Normalize FWHM to pixels
    beamFWHM /= gen_pixsize

    # Needed in this python version
    stdev = beamFWHM / (2*math.sqrt(2 * np.log(2)))

    print('beamfwhm = ', beamFWHM)
    print('stdev =',stdev)
    print('pixsize = ',pixsize)
    print('gen_pixsize = ',gen_pixsize)
    print('xgen = ',x_gen)
    print('ygen = ',y_gen)
    print('ioversamp =',ioversamp)
    print('ixcent = ', round(ixcent))
    print('iycent = ',round(iycent))
    print('npixx =', npixx)
    print('npixy =', npixy)
    # x_gen = int(x_gen - ixcent)
    # y_gen = int(y_gen - iycent)
    beamkern = Gaussian2DKernel(stdev,x_size = x_gen, y_size =  y_gen, mode = 'center')
    beamkern = np.array(beamkern)
    if ioversamp > 1:
        beamkern = 	rebin(beamkern,15,15)
    print('beamkern = ',beamkern)
# #   if we want this normalized than call with
#     if factor == 1:
# #       1D BEAM, this also now uses a scrpit called PSF_GAUSSIAN from astrolib
#         #convert FWHM to stdev\
#         stdev = beamFWHM / (2*math.sqrt(2*math.log(2)))
#         beamkern = scipy.signal.gaussian(x_gen, stdev) #this may work, but i have a feeling it won't work we might have to write our own function.
#     else:
#         pass #need to find a way to make a 2d gaussian
#
#         if ioversamp > 1:
#             # There doesnt seem to be a python equivalent
#             pass
#             #beamkern = rebin()

    plt.imshow(beamkern, interpolation='none', origin='lower')
    plt.xlabel('x [pixels]')
    plt.ylabel('y [pixels]')
    plt.colorbar()
    plt.show()
    print('beamkern = ', beamkern)
    return beamkern


def rebin(a, *args):
    shape = a.shape
    lenShape = len(shape)
    factor = np.asarray(shape)/np.asarray(args)
    evList = ['a.reshape('] + \
             ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
             [')'] + ['.mean(%d)'%(i+1) for i in range(lenShape)]
    print ''.join(evList)
    return eval(''.join(evList))


if __name__ == '__main__':
    get_spire_beam()
