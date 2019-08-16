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
from get_data import get_data
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
from gaussian import makeGaussian
from astropy.convolution import convolve_fft

def make_noise_mask(maps, col):
    """
    Inputs: Maps is the maps object that we get from get_data
            Col is either 0, 1, 2 to specify which band you want.
            0 = PSW, 1 = PMW, 2 = PLW
    """

    mapsize = maps[col]['error'].data.shape
    img = maps[col]['error'].data
    stddev = 24 / 2.355
    kernel = makeGaussian(8, 8, fwhm= 24)
    kernel = kernel / np.sum(kernel)
    # plt.imshow(kernel)
    # plt.show()

    kernel_full_size = np.zeros(mapsize)
    #find center of kernel
    xc = kernel_full_size.shape[0] / 2
    yc = kernel_full_size.shape[1] / 2
    x0 = xc - kernel.shape[0] / 2
    x1 = xc + kernel.shape[0] / 2
    y0 = yc - kernel.shape[1] / 2
    y1 = yc + kernel.shape[1] / 2
    kernel_full_size[x0:x1, y0:y1] = kernel
    # plt.imshow(kernel_full_size)
    # plt.show()

    fixed_image = convolve_fft(img, kernel_full_size)
    plt.imshow(fixed_image)
    plt.show()

    mask = np.empty(mapsize)
    for j in range(mapsize[0]):
        for k in range(mapsize[1]):
            if fixed_image[j,k] == 0:
                fixed_image[j,k] = 0.005
            if fixed_image[j,k] >= 0.004:
                mask[j,k] = 1
    return mask

if __name__ == '__main__':
    maps, err = get_data('a0370')
    mask = make_noise_mask(maps, 0)
    plt.imshow(mask)
    plt.show()
    # plt.imshow(mask)
    # plt.show()
