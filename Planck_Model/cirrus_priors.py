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
from clus_get_lambdas import *
from clus_get_data import *
from astropy.io import fits
from utilities import *
from make_map import *
import numpy as np
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
from astropy.wcs import WCS as world
import math
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve_fft as convolve


def one_map(name, map):
    #this is to call IRIS code
    ra_c = map['shead']['crval1']
    dec_c = map['shead']['crval12']
    coord = SkyCoord(ra=ra_c, dec=dec_c)
    b1950_coord = coord.transform_to(FK4(equinox='B1950'))
    ra_c_b = float(b1950_coord.ra.to_string(decimal=True))
    dec_c_b = float(b1950_coord.dec.to_string(decimal=True))

    name = map['clusname']
    #square patch of sky, this can be changed if desired.
    naxis = [map['shead']['naxis1'], map['shead']['naxis2']]

    head = make_header(map['pixsize'], naxis, ra_c_b, dec_c_b)

    iris_map = mosaic(head, band=4)

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

    iris_interp, r, d = interp_back_to_ref(iris_data, ra.ravel(), dec.ravel(), map['shead'])

    c       = 2.99792458e8            #m/s

    nu = c / clus_get_lambdas(map['band'], center=False)

    planck, ra, dec =  create_map(map['shead'], nu=nu)

    return planck, iris
