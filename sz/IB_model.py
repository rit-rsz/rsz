################################################################################
# NAME : IB_model.py
# DATE STARTED : Aug, 30, 2019
# AUTHORS : Dale Mercado & Benjamin Vaughan
# PURPOSE : Adding in the SZ isomap
# EXPLANATION : Create an iso-thermal beta model for testing.
# CALLING SEQUENCE :
# INPUTS : maps: Simmulated map objects
#          params: parameters from clus_get_params
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import numpy as np
from math import *
import os
import sys
sys.path.append('../utilities')
sys.path.append('../source_handling')
from config import *
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel
from astropy import units as u
import matplotlib.pyplot as plt

def IB_model(maps,params, verbose = 0):
    errmsg = False

    # Need to creat an array that will replicate the map size that we are using
    mapsize = maps['signal'].shape
    x = np.arange(0, mapsize[1])
    y = np.arange(0, mapsize[0])
    y = y[:, np.newaxis]
    # Need to find the center of the map
    ra = params['fidrad'] * u.deg
    dec = params['fidded'] * u.deg
    c = SkyCoord(ra, dec)
    # grabbing header from maps file
    # This is how it should be in the real compute rings
    w = wcs.WCS(maps['shead'])
    #converting ra/dec to pixel coords
    rc = params['rc']
    beta = params['beta']
    xcent, ycent = skycoord_to_pixel(c, w, origin=0)

    beta_map = (1. + (np.sqrt((x-xcent)**2 + (y-ycent)**2) / (rc / maps['pixsize']))**2)**((1. - 3. * beta) / 2.)
    return beta_map,errmsg
