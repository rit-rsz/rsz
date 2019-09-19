################################################################################
# NAME : make_noise_mask
# DATE STARTED : June 11, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE : A support routine for clus_compute_rings that makes a mask if noise is above a threshold.
#
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#           nsims (number of simulations/files to read in)
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../source_handling')
from clus_get_data import clus_get_data
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
from gaussian import makeGaussian
from astropy.convolution import convolve_fft
from gaussian import padGaussian

def clus_make_noise_mask(maps, col):
    """
    Inputs: Maps is the maps object that we get from clus_get_data
            Col is either 0, 1, 2 to specify which band you want.
            0 = PSW, 1 = PMW, 2 = PLW
    """

    mapsize = maps[col]['error'].shape
    img = maps[col]['error']
    stddev = 24 / 2.355
    kernel = makeGaussian(8, 8, fwhm= 4)
    kernel_full_size = padGaussian(img, kernel)

    #converting these to lists because the error message is annoying.
    img = img.tolist()
    kernel_full_size = kernel_full_size.tolist()
    fixed_image = convolve_fft(img, kernel_full_size)

    mask = np.zeros(mapsize, dtype=int)
    for j in range(mapsize[0]):
        for k in range(mapsize[1]):
            if fixed_image[j,k] == 0 or fixed_image[j,k] <= 1e-10 or np.isfinite(fixed_image[j,k]) == False:
                fixed_image[j,k] = 0.005
            if fixed_image[j,k] >= 0.004:
                mask[j,k] = int(1)

    ''' Temp fix for edge effects '''
    mask[0:10,:] = 1 # fix top 10 rows
    mask[-10:mapsize[0],:] = 1 # fix bottom 10 rows
    mask[:,0:10] = 1 # fix left 10 rows
    mask[:,-10:mapsize[1]] = 1 # fix right 10 rows
    return mask

if __name__ == '__main__':
    maps, err = clus_get_data('a0370')
    mask = make_noise_mask(maps, 0)
    # plt.imshow(mask,origin='lower')
    # plt.show()
