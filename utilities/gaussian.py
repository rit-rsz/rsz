################################################################################
# NAME : gaussian.py
# DATE STARTED : August 15, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE : Custom function for creating a gaussian, with option for a specifc location
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
#   Victoria Butler - 7/19 - edits to calfac for unit conversion
################################################################################
import numpy as np

def makeGaussian(x_size, y_size, fwhm = 3, center=None):
    """ Make a gaussian kernel.

    size is the length of a side of the square.
    fwhm is full-width-half-maximum.
    center is where you want center of gaussian located, default is center of array.
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
