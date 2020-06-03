################################################################################
# NAME : clus_interp_map.py
# DATE STARTED : May 20, 2020
# AUTHORS : Benjamin Vaughan
# PURPOSE : This script interpolates between SPIRE bands
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import numpy as np
import matplotlib.pyplot as plt
#from retrieve_iris_map import retrieve_iris_map
from scipy.interpolate import griddata
from astropy.wcs import WCS as world
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel

def interp_band_to_band(img, ref_map, map):

    map_size = map['signal'].shape
    ref_size = ref_map['signal'].shape

    #create a grid of RA/DEC coordinates for image we want to interpolate
    w = world(map['shead'])

    #holder for the x, y pixel coordinates that we want.
    x = np.arange(0, map_size[0])
    y = np.arange(0, map_size[1])
    xvec = np.repeat(x[:, np.newaxis], map_size[1], axis=1)
    yvec = np.repeat(y[np.newaxis, :], map_size[0], axis=0)

    c = pixel_to_skycoord(xvec, yvec, w, origin=0)
    #this is converting the pixel coords to right ascension and declination in fk4
    ra = np.asarray(c.ra.to_string(decimal=True), dtype=float)
    dec = np.asarray(c.dec.to_string(decimal=True), dtype=float)
    # reformat map data and coordinates
    data = np.ravel(img)
    points = np.column_stack((np.ravel(ra), np.ravel(dec)))

    #create a agrid of RA/DEC coordinates that we want to interpolate over
    ref_w = world(ref_map['shead'])
    ref_grid_x, ref_grid_y = np.mgrid[0:ref_size[0], 0:ref_size[1]]
    ref_grid_ra, ref_grid_dec = ref_w.wcs_pix2world(ref_grid_x, ref_grid_y, 0)

    #do the interpolation
    interp_map = griddata(points, data, (ref_grid_ra, ref_grid_dec), method='linear')
    return interp_map
