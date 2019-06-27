################################################################################
# NAME : starfindertest.py
# DATE STARTED : June 26, 2019
# AUTHORS : Dale Mercado
# PURPOSE : This is just used to test the different modules that could potentially
#           replace starfinder.
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import matplotlib as plt
import matplotlib.pyplot as p
import sys
sys.path.append('source_handling')
sys.path.append('utilities')
from get_data import *
import numpy as np
from writefits import *
import config

def starfindertest(clusname):
    # make a table of Gaussian sources
    from astropy.table import Table
    table = Table()
    table['amplitude'] = [50, 70, 150, 210]
    table['x_mean'] = [160, 25, 150, 90]
    table['y_mean'] = [70, 40, 25, 60]
    table['x_stddev'] = [15.2, 5.1, 3., 8.1]
    table['y_stddev'] = [2.6, 2.5, 3., 4.7]
    table['theta'] = np.array([145., 20., 0., 60.]) * np.pi / 180.

    # make an image of the sources with Gaussian noise
    from photutils.datasets import make_gaussian_sources_image
    from photutils.datasets import make_noise_image
    shape = (100, 200)
    sources = make_gaussian_sources_image(shape, table)
    noise = make_noise_image(shape, type='gaussian', mean=0.,
                             stddev=5., random_state=12345)
    image = sources + noise

    # detect the sources
    from photutils import detect_threshold, detect_sources
    threshold = detect_threshold(image, snr=3)
    from astropy.convolution import Gaussian2DKernel
    sigma = 3.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))   # FWHM = 3
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    segm = detect_sources(image, threshold, npixels=5,
                          filter_kernel=kernel)

    # plot the image and the segmentation image
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
    ax1.imshow(image, origin='lower', interpolation='nearest')
    ax2.imshow(segm.data, origin='lower', interpolation='nearest')

if __name__ == '__main__':
    starfindertest('rxj1347')