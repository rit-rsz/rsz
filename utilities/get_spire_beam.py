
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
import numpy as np
import matplotlib.pyplot as plt
from math import *
from astropy.io import fits
import scipy.signal
import os
from astropy.convolution import Gaussian2DKernel, Gaussian1DKernel
from get_spire_beam_fwhm import *



def get_spire_beam(band=None, pixsize=0,npixx=0, npixy=0,
                   xcent=0, ycent=0,bolometer=0, fwhm='',
                   norm=0, oversamp=0, verbose=1,
                   factor=0):
    errmsg = False
    if bolometer != 0:
        if verbose:
            print('PSFs for specific bolometer not yet supported -- ignoring')

    # Check if we have been given a band
    # If not make an assumption
    if band == None:
        if verbose:
            print('Band parameter not supplied, assuming PSW')
        band = 'PSW'
        # Check if we have been givin a pixel size
    if pixsize == 0:
        # band = upper(band)
        # units arcsec/pixel
        if band == 'PSW':
            pixsize = 6
        if band == 'PMW':
            pixsize = 8. + (1/3)
        if band == 'PLW':
            pixsize = 12
        else:
            print('Unknown band'+band)
        if verbose:
            print('pixsize paramter not supplied assuming %s arcsec' , (pixsize))


    if len(fwhm) == 0:
        beamFWHM = get_spire_beam_fwhm(band)
    else:
        beamFWHM = fwhm

    if beamFWHM < 0:
        print('Invalid Beam FWHM value'+ str(beamFWHM))

    # Check if we've been given the map size, if not assume something
    # npixx/npixy will be the final number of pixels
    if npixx == 0:
        npixx = round(beamFWHM * 5 / (pixsize))
        if npixx % 2 != 1:
            npixx +=1
            if verbose:
                print('npixx not supplied, using ",I0"')
    #If no y size then assume same as x
    if npixy == 0:
        npixy = npixx

    #Make sure that these have been cast from a float properly or we get errors
    npixx = int(ceil(npixx))
    npixy = int(ceil(npixy))

    if npixx % 2 != 1 and verbose:
        print('WARNING: npixx not odd, so PSF will not be centered')
    if npixy % 2 != 1 and verbose:
        print('WARNING: npixy not odd, so psf will not be centered')

    # Now deal with oversampling
    if oversamp == 0:
        ioversamp = 7
    else:
        ioversamp = round(oversamp) #in case user provides float
        if ioversamp % 2 != 1:
            print('Oversamp must be an odd interger!')

    x_gen = npixx * ioversamp
    y_gen = npixy * ioversamp
    gen_pixsize = np.float64(pixsize) / ioversamp

    # Check if we have been givin the center, if not assume middle
    if xcent == 0:
        # if verbose:
        #     print('xcent parameter not supplied, assuming array center')
        ixcent = x_gen / 2
    else:
        # Adjust for oversampling
        ixcent = xcent * ioversamp
        if ioversamp > 1:
            ixcent = ixcent + ioversamp / 2

    if ycent == 0:
        # if verbose:
        #     print('ycent parameter not supplied, assuming array center')
        iycent = y_gen / 2
    else:
        iycent = ycent * ioversamp
        if ioversamp > 1:
            iycent = iycent +ioversamp / 2

    # Normalize FWHM to pixels
    beamFWHM /= gen_pixsize
    # Convert the FWHM to a standard deviation for astropy fit. From psf_gaussian
    stdev = beamFWHM / (sqrt(8 * log(2)))

    # If we want this normalized then call with norm flag set
    if factor:
        # 1D beam
        beamkernraw = Gaussian1DKernel(stdev, x_size=x_gen)
        beamkern = np.array(beamkernraw)
        if norm:
            beamkern = beamkern / beamkern.max()
        if ioversamp > 1:
            beamkern = rebin(beamkern,(npixx,npixy))
    else:
        beamkernraw = Gaussian2DKernel(stdev,x_size = x_gen, y_size =  y_gen)
        beamkern = np.array(beamkernraw)
        if norm:
            beamkern = beamkern / beamkern.max()
        if ioversamp > 1:
            beamkern = 	rebin(beamkern,(npixx,npixy))



    # # Use for debugging
    # plt.plot for 1d
    # plt.plot(beamkern, drawstyle='steps')
    # plt.imshow for 2d & colorbar
    # plt.imshow(beamkern, interpolation='none', origin='lower')
    # plt.colorbar()
    # plt.xlabel('x [pixels]')
    # plt.ylabel('y [pixels]')
    # plt.show()
    # beamkern = 1

    return beamkern


def rebin(a, new_shape):
    shape = a.shape
    M = int(shape[0])
    N = int(shape[1])
    m, n = new_shape
    if m<M:
        return a.reshape((m,int(M/m),n,int(N/n))).mean(3).mean(1)
    else:
        return np.repeat(np.repeat(a, m/M, axis=0), n/N, axis=1)

if __name__ == '__main__':
    get_spire_beam(norm=1,factor = 0)
