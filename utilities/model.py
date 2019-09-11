# NAME : xid_model.py
# DATE STARTED : July 19, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE : a script to create a mask map from the XID info.
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import json
import sys
import matplotlib.pyplot as plt
from astropy.io import fits
sys.path.append('../source_handling')
print(sys.path)
import numpy as np
import os
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.modeling.functional_models import Gaussian2D
from scipy.signal import convolve as convolver
from photutils.psf import IntegratedGaussianPRF
from math import *
from photutils import CircularAperture
from astropy.wcs.utils import pixel_to_skycoord
from astropy.stats import sigma_clipped_stats
from photutils import datasets
from photutils import DAOStarFinder

def create_model(maps, cats):
        bands = [18, 25, 36]
        gal_clusts = []
        for i in range(len(cats[i]['sflux'])):
            print('generating mask for %s at %s' % (maps[i]['name'], maps[i]['band']))
            fluxes = cats[i]['sflux']
            naxis = maps[i]['astr']['NAXIS']
            gal_clust = np.zeros( naxis)
            hdul = fits.open(maps[i]['file'])
            w = WCS(hdul[1].header)
            ra = np.array(cats[i]['sra']) * u.deg
            dec = np.array(cats[i]['sdec']) * u.deg
            c = SkyCoord(ra, dec)
            px, py = skycoord_to_pixel(c, w, 1)
            y_size = hdul[1].data.shape[0]
            x_size = hdul[1].data.shape[1]
            for j in range(len(cats[i][''])):
                kern = makeGaussian(x_size=x_size, y_size=y_size, fwhm =bands[i] / maps[i]['pixsize'], center=(px[j],py[j]))
                kern = np.asarray(kern)
                kern = kern / kern.max()
                coefficient = fluxes[j]
                psf = kern * coefficient
                gal_clust = gal_clust + psf
            gal_clusts.append(gal_clust)
            print('finished generating mask for %s' % (maps[i]['band']))
        return gal_clusts

def makeGaussian(x_size, y_size, fwhm = 3, center=None):
    """ Make a square gaussian kernel.

    size is the length of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """

    sigma = fwhm / 2.355
    x = np.arange(0, x_size, 1, float)
    y = np.arange(0, y_size, 1, float)
    y = y[:,np.newaxis]

    if center is None:
        x0 = x_size // 2
        y0 = y_size // 2
    else:
        x0 = center[0]
        y0 = center[1]

    sigma = fwhm / 2.355

    return np.exp(-1 * ((x-x0)**2 + (y-y0)**2) / (2*sigma**2))
