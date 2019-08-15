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

def make_noise_mask(maps, col):
    """
    Inputs: Maps is the maps object that we get from get_data
            Col is either 0, 1, 2 to specify which band you want.
            0 = PSW, 1 = PMW, 2 = PLW
    """

    mapsize = maps[col]['error'].data.shape
    img = maps[col]['error'].data
    #
    # plt.imshow(img, vmin=0.002, vmax=0.005)
    # plt.show()
    plt.imshow(maps[col]['mask'])
    plt.show()
    stddev = 24 / 2.355
    kernel = Gaussian2DKernel(x_stddev=stddev)
    fixed_image = interpolate_replace_nans(img, kernel)
    # plt.imshow(fixed_image)
    # plt.show()
    mask = np.empty(mapsize)
    for j in range(mapsize[0]):
        for k in range(mapsize[1]):
            if fixed_image[j,k] == 0:
                fixed_image[j,k] = 0.005
    #call to make a gaussian function and convol it with the map to smooth it out, but I don't think we need that?
            if fixed_image[j,k] >= 0.004:
                mask[j,k] = 1
    return mask

if __name__ == '__main__':
    maps, err = get_data('a2218')
    mask = make_noise_mask(maps, 0)
    plt.imshow(mask)
    plt.show()
