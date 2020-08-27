################################################################################
# NAME : make_map.py
# DATE STARTED : Julu 14, 2020
# AUTHORS : Benjamin Vaughan
# PURPOSE : To create a data structure to implement into PCAT
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import sys
sys.path.append('../rsz/utilities')
sys.path.append('../iris/')
from utility import *
from mosaic import *
from clus_get_lambdas import *
from astropy.io import fits
from astropy.coordinates import FK4
from make_map import *
import numpy as np
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
from astropy.wcs import WCS as world
import math
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve_fft as convolve
from planck_utilities import *


def one_map(name, header, band):

    if 'CDELT1' in header.keys():
        #calculating the pixsize based off of astrometry parameters.
        pixsize = 3600 * np.mean([abs(header['CDELT1']),abs(header['CDELT2'])])
        #if the given astrometry has CDELT values and not a cd matrix create a cd matrix.
        header['cd_1'] = header['CDELT1']
        header['cd_2'] = 0
        header['cd_1'] = 0
        header['cd_2'] = header['CDELT2']


    else:
        #if the astrometry is given with a cd matrix calculate the pixsize this way.
        pixsize = 3600 * \
                  np.mean([abs(header['CD1_1']+  header['CD2_1']), \
                        abs(header['CD2_1'] + header['CD2_2'])])


    #this is to call IRIS code
    ra_c = header['crval1'] * u.deg
    dec_c = header['crval2'] * u.deg
    cent = [ra_c.value, dec_c.value]
    coord = SkyCoord(ra=ra_c, dec=dec_c)
    b1950_coord = coord.transform_to(FK4(equinox='B1950'))
    ra_c_b = float(b1950_coord.ra.to_string(decimal=True))
    dec_c_b = float(b1950_coord.dec.to_string(decimal=True))
    #square patch of sky, this can be changed if desired.
    naxis = [header['naxis1'], header['naxis2']]

    iris_head = make_header(pixsize, naxis, ra_c_b, dec_c_b)
    iris_map = mosaic(iris_head, band=4)

    x = np.arange(0, iris_head['NAXIS1'])
    y = np.arange(0, iris_head['NAXIS2'])
    yvec = np.repeat(x[:, np.newaxis], iris_head['NAXIS2'], axis=1)
    xvec = np.repeat(y[np.newaxis, :], iris_head['NAXIS1'], axis=0)
    w = world(iris_head)
    c = pixel_to_skycoord(xvec, yvec, w, origin=0)
    c = c.transform_to('icrs')

    #this is converting the pixel coords to right ascension and declination in fk4
    ra = np.asarray(c.ra.to_string(decimal=True), dtype=float)
    dec = np.asarray(c.dec.to_string(decimal=True), dtype=float)

    iris_interp, r, d = interp_back_to_ref(iris_map, ra.ravel(), dec.ravel(), header)

    c       = 2.99792458e8            #m/s

    nu = c / clus_get_lambdas(band, center=False)

    planck, ra, dec, avg_T, avg_beta, avg_tau =  create_map(header, nu=nu)
    planck_head = save_header(header, avg_beta, avg_T, avg_tau, 'Planck_' + band, 'Jy/Beam', cent)
    return planck, iris_interp, avg_T, avg_beta, avg_tau, planck_head
