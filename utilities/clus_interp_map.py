################################################################################
# NAME : clus_interp_map.py
# DATE STARTED : May 20, 2020
# AUTHORS : Benjamin Vaughan
# PURPOSE : This script interpolates IRAS coordinates to SPIRE coordinates
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


def interp_map(maps, params):

    #Teresa's script seems to use a corner and not a center so might need to change this
    ra_cent = params['nomrad']
    dec_cent = params['nomded']

    # Create IRIS map at same coordinates
    retrieve_iris_map(ra,dec)
    irismap = fits.open('IRIS/iris_fx.fits')
    irisra = fits.open('IRIS/iris_ra.fits')
    irisdc = fits.open('IRIS/iris_dc.fits')

    # Get IRIS data and coordinates
    iris_data = np.ravel(irismap[0].data)
    iris_points = np.column_stack((np.ravel(irisra[0].data), np.ravel(irisdc[0].data)))

    # Close IRIS files
    irismap.close()
    irisra.close()
    irisdc.close()

    for i in range(len(maps)):
        signal = maps[i]['signal']
        w = wcs.WCS()
        #create spire grid
        map_size = signal.shape
        grid_x, grid_y = np.mgrid[0:map_size[0],0:map_size[1]]
        grid_ra, grid_dec = w.wcs_pix2world(grid_x, grid_y, 0)

        # Create interpolated IRIS image
        iris_interp = griddata(iris_points, iris_data, (grid_ra, grid_dec), method='linear')
        return iris_interp

def interp_band_to_band(ref_map, map):

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
    data = np.ravel(map['signal'])
    points = np.column_stack((np.ravel(ra), np.ravel(dec)))

    print(data.shape)
    print(map['signal'].shape)

    #create a agrid of RA/DEC coordinates that we want to interpolate over
    ref_w = world(ref_map['shead'])
    ref_grid_x, ref_grid_y = np.mgrid[0:ref_size[0], 0:ref_size[1]]
    ref_grid_ra, ref_grid_dec = ref_w.wcs_pix2world(ref_grid_x, ref_grid_y, 0)

    #do the interpolation
    interp_map = griddata(points, data, (ref_grid_ra, ref_grid_dec), method='linear')
    print(interp_map.shape)
    return interp_map
