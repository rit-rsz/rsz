################################################################################
# NAME : SIDES_TEST
# DATE STARTED : January 17, 2020
# AUTHORS : Benjamin Vaughan
# PURPOSE : As part of the development of the blank field testing with the SIDES
# catalog we need a model for picking out squares of a similar size to our SPIRE
# maps from the 2 deg^2 SIDES map that don't have significant overlap.
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
import time
from math import *
from astropy.io import fits
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
from astropy.wcs import WCS as world
from astropy.coordinates import SkyCoord
import astropy.units as u
import sys
sys.path.append('../source_handling')
sys.path.append('../utilities')
from clus_get_data import *
from read_sides import ret_sides

def create_catalogues(pixsize, cutout, sides_size, nsamples):
    '''
    Pixsize : the pixel size
    cutout  : tuple, array, or list : the length of the sides of the cutout
    sides_size : area of sides map in square degrees (should be 2)
    nsamples : the number of cutouts that we want to maek
    '''
    x_size = y_size = sqrt(sides_size) #calculate lengths from the total area
    p_x_size = x_size * 3600 / pixsize
    p_y_size = y_size * 3600 / pixsize #converting the sides catalogue to pixel scales
    sides_p_size = [p_x_size, p_y_size] #sides catalogue in pixels

    #find the central values for our cut outs.
    x, y = map_selector(sides_p_size, cutout, nsamples)

    c_pix = [[xi, yi] for xi, yi in zip(x,y)]
    print(c_pix, 'center pixel')
    return c_pix

def map_selector(m_size, s_size, n_samp):
    """
    m_size : catalog size
    s_size : cut out size
    n_samp : number of cuttouts
    plot   : plotting flag, True means plot and False means don't plot
    """
    #create a new map that allows us to select center pixels for potential SPIRE
    #maps that don't overextend the edges of the SIDES map.

    #calculate needed floating point values to generate centers within the correct
    #area
    half_f = [floor(s_size[0] / 2), floor(s_size[1] / 2)]
    n_shape = [m_size[0] - half_f[0], m_size[1] - half_f[1]]

    #create a square grid of our central values
    num_x = num_y = int(sqrt(n_samp))
    x = np.linspace(half_f[0], n_shape[0], num_x, dtype=np.float32)
    y = np.linspace(half_f[1], n_shape[1], num_y, dtype=np.float32)
    coords = np.stack(np.meshgrid(x, y), -1).reshape(-1,2)
    px = coords[:,0]
    py = coords[:,1]
    return px, py

def populate_cutouts(sides_catalogue, c_pix, pixsize, band, cutout):
    if band == 'PSW':
        b = 3
    elif band == 'PMW':
        b = 4
    elif band == 'PLW':
        b = 5

    py = sides_catalogue[1] #y values in the sides catalogue
    px = sides_catalogue[0] #x values in the sides catalogue
    f  = sides_catalogue[b] #flux values at our specified band
    z  = sides_catalogue[2] #red shift in the sides catalogue

    #find splices for the cutout
    cx = c_pix[0] * pixsize / 3600.
    cy = c_pix[1] * pixsize / 3600.
    sx_2 = cutout[0] / 2. * pixsize / 3600.
    sy_2 = cutout[1] / 2. * pixsize / 3600.

    ux = cx + sx_2
    dx = cx - sx_2
    uy = cy + sy_2
    dy = cy - sy_2

    print(ux, dx, uy, dy)
    #cutouts of map
    a = time.time()
    good_x = np.where(np.logical_and(px >= dx, px <= ux))[0]
    good_y = np.where(np.logical_and(py >= dy, py <= uy))[0]
    good_xy = np.intersect1d(good_x, good_y)
    good_f = np.where(f[good_xy] > 0)
    good = np.intersect1d(good_xy, good_f)
    b = time.time()
    print('Number of sources found: %s in %s' % (len(f[good]), b-a))
    return px[good], py[good], f[good], z[good]


if __name__ == '__main__':
    # map_selector((50,50), (5,5), 100)
    a = time.time()
    maps = clus_get_data('rxj1347', 1, manpath=0, resolution = 'nr', bolocam=None,
                        verbose = 1, version = '1', manidentifier=None, sgen=None, nsim=0, testflag=0)[0]
    bands = ['PSW', 'PMW', 'PLW']
    master_list = ret_sides()
    for i in range(len(maps)):
        pixsize = maps[i]['pixsize']
        # cutout = maps[i]['signal'].shape
        cutout = [290, 270]
        band = maps[i]['band']
        c_pix_list = create_catalogues(pixsize, cutout, 2, 100)
        j = 0
        for c in c_pix_list:
            filename = '../sides_sims/sides_%s_sim%s' % (band, j)
            x, y, f, z = populate_cutouts(master_list, c, pixsize, band, cutout)
            truth_table = {'RA' : x.to_numpy('float'),
                           'DEC': y.to_numpy('float'),
                           'Flux' : f.to_numpy('float'),
                           'Redshift' : z.to_numpy('float')}
            # np.save(filename, truth_table, allow_pickle=True)
            j += 1
    # b = time.time()
    # c_pix_list = create_catalogues(1, [2*3600,2*3600], 2, 1)
    # for j in range(len(maps)):
    #     for c in c_pix_list:
    #         filename = '../sides_sims/SIDES_FULL_FIELD_%s' % (bands[j])
    #         x, y, f, z = populate_cutouts(master_list, c, 1, bands[j], [2 * 3600, 2 * 3600])
    #         truth_table = {'RA' : x,
    #                        'DEC': y,
    #                        'Flux' : f,
    #                        'Redshift' : z}
    #         np.save(filename, truth_table, allow_pickle=True)

    print('Total Runtime :', b-a, "Start :", a, "Finish :", b)
    # plt.show()


    # populate_cutouts_test((50,50), 5000, (10,15), (5,5))