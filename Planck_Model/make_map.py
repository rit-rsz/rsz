################################################################################
# NAME : make_map.py
# DATE STARTED : June 25, 2020
# AUTHORS : Benjamin Vaughan
# PURPOSE : the purpose of this program is to make a map of Cirrus Emission using
# the results of the Dust Fitting Analysis done by the Planck Team on their 2015
# GNILC-Dust models.
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
from astropy.io import fits
import numpy as np
from math_functions import *
from scipy.interpolate import griddata
from astropy.wcs import WCS as world
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
import matplotlib.pyplot as plt
from astropy_healpix import HEALPix
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy.interpolate import griddata
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
import planck_config

def create_map(ref_head, nu):
    '''
    Inputs : filenames (str array) - the paths to and the filenames of the component maps in order of Tau, Temp, Beta (see example.py)
             ref_head  (astropy.header) - the header for the reference map to interpolate onto
             ref_pixsize (float) - the pixel size of the reference map
             ref_mapsize (int array) - a set of the values for the length of each axes of the map
             center (float array) - a set of the values for the center of the field of interest
    Outputs : Interped_map (float array) - an array of fluxes interpolated onto the reference map's grid.
              ragrd (float array) - an array of the Right Ascension values used in the interpolation
              decgrd (float array) - an array of the Declination values used in the interpolation
    Notes : The Planck model has an offset in intensity due to the fact that in Planck Collaborationet al. 2016 correlated their
            fit to regions where H1 is present and therefore had an offset in their fit. For more details see their paper.
    '''

    tau_name = config.PLANCKDATA + 'COM_CompMap_Dust-GNILC-Model-Opacity_2048_R2.01.fits'
    temp_name = config.PLANCKDATA + 'COM_CompMap_Dust-GNILC-Model-Temperature_2048_R2.01.fits'
    beta_name = config.PLANCKDATA + 'COM_CompMap_Dust-GNILC-Model-Spectral-Index_2048_R2.01.fits'

    filenames = [tau_name, temp_name, beta_name]


    if 'CD1_1' in ref_head.keys():
        ref_pixsize = 3600 *  np.mean([abs(ref_head['CD1_1'] + ref_head['CD2_1']), abs(ref_head['CD2_1'] + ref_head['CD2_2'])])
    elif 'CDELT1' in ref_head.keys():
        ref_head['CD1_1'] = - 1 * ref_head['CDELT1']
        ref_head['CD1_2'] = 0
        ref_head['CD2_1'] = 0
        ref_head['CD2_2'] = ref_head['CDELT2']
        ref_pixsize = 3600 *  np.mean([abs(ref_head['CD1_1'] + ref_head['CD2_1']), abs(ref_head['CD2_1'] + ref_head['CD2_2'])])
    ref_mapsize = [ref_head['NAXIS2'], ref_head['NAXIS1']]
    center = [ref_head['CRVAL1'], ref_head['CRVAL2']]

    #this is to oversize the map so that the data has a larger area than the interpolating grid to reduce edge effects
    ref_head['naxis1'] = ref_head['naxis1'] * 2
    ref_head['naxis2'] = ref_head['naxis2'] * 2
    ref_head['crpix1'] = ref_head['crpix1'] * 2
    ref_head['crpix2'] = ref_head['crpix2'] * 2


    I_to_MJy = 1e20 #converts from standard units of Intensity to MJy
    param_values = [] #Beta = 2, Temperature = 1, Tau = 0
    for f in filenames:
        data, pixsize, x_side, y_side, ra, dec = read_in_fits(f, center, ref_head, ref_pixsize, ref_mapsize)
        param_values.append(data)
    I = calc_intensity(param_values[2], param_values[0], param_values[1], nu)

    I_map = np.reshape(I, (x_side, y_side)) * I_to_MJy


    #this fixes the astrometry in the map to account for the oversizing done previously
    ref_head['naxis1'] = ref_head['naxis1'] / 2
    ref_head['naxis2'] = ref_head['naxis2'] / 2
    ref_head['crpix1'] = ref_head['crpix1'] / 2
    ref_head['crpix2'] = ref_head['crpix2'] / 2

    interped_map, ragrd, decgrd = interp_back_to_ref(I_map, ra, dec, ref_head)


    return interped_map, ragrd, decgrd
    # return interped_map

def read_in_fits(filename, center, ref_head, ref_pixsize=8, ref_mapsize=260):
    '''
    Purpose : This function reads in the fits files for the components and parses them so that we are left with data only for our field.
    Inputs: filename (str) - the singular filename for a component used in the MBB fit
            center (float array) - center of the field of interest
            ref_head (Astropy.header) - header for the reference field
            ref_pixsize - pixel size of the reference map
            ref_mapsize - mapsize of the reference map
    Outputs: map (float array) - an array of flux values at the given field
             pixsize (float) - pixel size of the uninterpolated component maps
             x_side (int) - the length of the x-axis of the map
             y_side (int) - the length of the y-axis of the map
             RA_grid (float array) - grid of the Right Ascension values used to pull out components
             DEC_grid (float array) - grid of the Declination values used to pull out components
    '''

    hdul = fits.open(filename)
    head = hdul[1].header
    if 'Temperature' in filename:
        data = hdul[1].data.field('TEMP')
        error = hdul[1].data.field('ERR_TEMP')
    elif 'Spectral-Index' in filename:
        data = hdul[1].data.field('BETA')
        error = hdul[1].data.field('ERR_BETA')

    elif 'Opacity' in filename:
        data = hdul[1].data.field('TAU353')
        error = hdul[1].data.field('ERR_TAU')

    else:
        data = hdul[1].data.field(0)
    nside = head['NSIDE']
    order = head['ORDERING']
    hdul.close()

    #Galactic Coordinate System
    hp = HEALPix(nside=nside, order=order, frame='galactic')
    #create a pixel grid in terms of the nu=353 grid for GNILC to create our intensity maps
    pixsize = hp.pixel_resolution.to(u.arcsecond).value

    #* 10 is for boosting the size of the map so to fix edge effects from interpolation
    map_arc_x = ref_mapsize[0] * 2 * ref_pixsize #map size in arcseconds
    map_arc_y = ref_mapsize[1] * 2 * ref_pixsize

    npixxside = ceil(map_arc_x / pixsize) #convert to map size in pixels for nu = 353 map.
    npixyside = ceil(map_arc_y / pixsize)

    #* 10 is for boosting the size of the map so to fix edge effects from interpolation
    x  = np.linspace(0, ref_mapsize[0] * 2,   npixxside)
    y  = np.linspace(0, ref_mapsize[1] * 2,   npixyside)

    X, Y = np.meshgrid(x, y)
    w = world(ref_head)
    skycoords = pixel_to_skycoord(X.ravel(), Y.ravel(), wcs=w, origin=0)
    RA_grid = np.asarray(skycoords.ra.to_string(decimal=True), dtype='float') * u.deg
    DEC_grid = np.asarray(skycoords.dec.to_string(decimal=True), dtype='float') * u.deg




    # coords = SkyCoord(RA_grid.ravel(), DEC_grid.ravel(), frame='icrs')
    coords = SkyCoord(ra=RA_grid.ravel(), dec=DEC_grid.ravel(), frame='icrs')
    gal_coords = coords.galactic

    map = hp.interpolate_bilinear_skycoord(gal_coords, data)

    x_side = len(x)
    y_side = len(y)
    return map, pixsize, y_side, x_side, RA_grid, DEC_grid


def interp_back_to_ref(img, ra, dec, ref_head):
    '''
    Purpose: Perform the inperolation of the completed Planck Map to a reference header
    Inputs: img - the Planck flux map
            ra  - grid of the Right Ascension values used to create the Planck Map
            dec - grid of the Decliation values used to create the Planck Map
            ref_head - the reference header for the interpolatino grid
            ref_shape - the shape of the interpolation grid
    '''

    ref_shape = [ref_head['naxis1'], ref_head['naxis2']]

    map_size = img.shape

    # reformat map data and coordinates
    data = np.ravel(img)
    points = np.column_stack((np.ravel(ra), np.ravel(dec)))

    #create a agrid of RA/DEC coordinates that we want to interpolate over
    ref_w = world(ref_head)
    ref_grid_x, ref_grid_y = np.mgrid[0:ref_shape[0], 0:ref_shape[1]]
    ref_grid_ra, ref_grid_dec = ref_w.wcs_pix2world(ref_grid_x, ref_grid_y, 0)

    #do the interpolation
    interp_map = griddata(points, data, (ref_grid_ra, ref_grid_dec), method='linear')
    final_map = np.swapaxes(interp_map, 0, 1)
    return final_map, ref_grid_ra, ref_grid_dec
