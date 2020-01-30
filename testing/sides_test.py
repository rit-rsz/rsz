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
from math import *
from astropy.io import fits
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
from astropy.wcs import WCS as world
from astropy.coordinates import SkyCoord
import astropy.units as u
import sys
sys.path.append('../source_handling')
from clus_get_data import *

def create_catalogues(maps, sides_size):
    x_size = y_size = sqrt(sides_size)
    for i in range(len(maps)):
        pixsize = maps[0][i]['pixsize']
        p_x_size = x_size * 3600 / pixsize
        p_y_size = y_size * 3600 / pixsize
        cutout_size = maps[0][i]['signal'].shape
        sides_p_size = [p_x_size, p_y_size]

        x, y = map_selector(sides_p_size, cutout_size, 100)
        c_pix = np.asarray(zip(x,y))


def map_selector(m_size, s_size, n_samp):
    """
    m_size : catalog size
    s_size : cut out size
    n_samp : number of cuttouts
    """
    print(s_size, m_size)
    #create a new map that allows us to select center pixels for potential SPIRE
    #maps that don't overextend the edges of the SIDES map.

    #calculate needed floating point values to generate centers within the correct
    #area
    half_f = [floor(s_size[0] / 2), floor(s_size[1] / 2)]
    n_shape = [m_size[0] - half_f[0], m_size[1] - half_f[1]]
    #generate random center pixels with x, y coordinates
    # c_pix = []
    # for i in range(n_samp):
    #     cx = np.random.random() * n_shape[0] + floor(s_size[0] / 2)
    #     cy = np.random.random() * n_shape[1] + floor(s_size[1] / 2)
    #     c_pix.append([cx, cy])
    # # final_ds(c_pix)
    # c_pix = np.asarray(c_pix)

    #for generating a grid instead
    #create a square grid
    num_x = num_y = int(sqrt(n_samp))
    x = np.linspace(half_f[0], n_shape[0], num_x, dtype=np.float32)
    y = np.linspace(half_f[1], n_shape[1], num_y, dtype=np.float32)
    coords = np.stack(np.meshgrid(x, y), -1).reshape(-1,2)
    plt.scatter(coords[:,0], coords[:,1], label='coordinate center')
    for i in range(n_samp):
        c_pix = coords[i,0], coords[i,1]
        plot_squares(c_pix, s_size)
    #plot the greater map
    g_ux = m_size[0] #upper x
    g_dx = 0 #lower x
    g_uy = m_size[1] #upper y
    g_dy = 0 #lower y
    #coordinates for grid
    gx = [g_dx, g_dx, g_ux, g_ux, g_dx]
    gy = [g_dy, g_uy, g_uy, g_dy, g_dy]

    #show that the centers are within an area such that edges of the map don't
    #overextend the SIDES catalog
    #find the inner corners
    i_ux = 0 + s_size[0] / 2.
    i_uy = 0 + s_size[1] / 2.
    i_dx = m_size[0] - s_size[0] / 2.
    i_dy = m_size[1] - s_size[1] / 2.

    ix = [i_dx, i_dx, i_ux, i_ux, i_dx]
    iy = [i_dy, i_uy, i_uy, i_dy, i_dy]


    plt.plot(ix, iy, c='green', label='inner square')

    plt.plot(gx, gy, c='red', label='outer square')
    plt.show()
    px = coords[:,0]
    py = coords[:,1]
    return px, py

def plot_squares(c_pix, s_size):
    """
    c_pix : center pix coordinates
    s_size: cut out size
    """
    #Test code for overlapping maps
    #break up inputs into individual values
    cx = c_pix[0]
    cy = c_pix[1]
    sx_2 = s_size[0] / 2.
    sy_2 = s_size[1] / 2.
    #find the corner pieces of the square
    ux = cx + sx_2
    dx = cx - sx_2
    uy = cy + sy_2
    dy = cy - sy_2
    #create a list of points to draw the grid
    ix = [dx, dx, ux, ux, dx]
    iy = [dy, uy, uy, dy, dy]
    #plot the grid
    plt.plot(ix, iy, c='blue')

def populate_cutouts(sides_catalogue, nsamples, c_pix, s_size):
    #####code for testing
    # x = np.random.random_sample(nsamples)
    # y = np.random.random_sample(nsamples)
    # px = [xi * m_shape[0] for xi in x]
    # py = [yi * m_shape[1] for yi in y]
    # px = np.asarray(px)
    # py = np.asarray(py)
    # plt.scatter(px, py)
    # plt.title('Test Source Field')
    # plt.show()
    #find splices for the cutout
    cx = c_pix[0]
    cy = c_pix[1]
    sx_2 = s_size[0] / 2.
    sy_2 = s_size[1] / 2.
    ux = cx + sx_2
    dx = cx - sx_2
    uy = cy + sy_2
    dy = cy - sy_2

    #splice of the map would be: cx-sx_2:cx+sx_2, cy-sy_2:cy+sy_2
    #cutouts of map
    good_x = np.where(np.logical_and(px >= dx, px <= ux))[0]
    good_y = np.where(np.logical_and(py >= dy, py <= uy))[0]
    good = [gx for gx in good_x if gx in good_y]
    return px[good], py[good]



##########################################Old version of generating cutouts might want this later
# def calc_ds(c_list, c_pix):
#     #calculate manhattan distance and return a list of distances
#     d_list = []
#     for c in c_list:
#         d = abs(c[0] - c_pix[0]) + abs(c[1] - c_pix[1])
#         if d < 2:
#             d_list.append(d)
#     return d_list
#
# def final_ds(c_list):
#     #check the manhattan distances for all points
#     ds = []
#     i = 0
#     for c1 in c_list:
#         i += 1
#         j = 0
#         for c2 in c_list:
#             j += 1
#             d = abs(c1[0] - c2[0]) + abs(c1[1] - c2[1])
#             print('object %s with coordinates %s compared to object %s with coordinates %s' % (i, c1, j, c2), d)
#
# def check_c(c_list, c_pix, n_shape, s_size):
#     iter = 0
#     if len(c_list) > 0:
#         d_list = calc_ds(c_list, c_pix)
#         while len(d_list) > 0:
#             if iter > 100:
#                 print('Error: exceeded iteration limit.')
#                 break #break the while loop if it goes over 100 iterations
#             iter += 1
#             cx = np.random.random() * n_shape[0] + floor(s_size[0] / 2)
#             cy = np.random.random() * n_shape[1] + floor(s_size[1] / 2)
#             d_list = calc_ds(c_list, [cx, cy])
#         print('check')
#         return c_pix
#     else:
#         return c_pix

if __name__ == '__main__':
    # map_selector((50,50), (5,5), 100)
    maps = clus_get_data('rxj1347', manpath=0, resolution = 'nr', bolocam=None,
                verbose = 1, version = '1', manidentifier=None, sgen=None, nsim=0, testflag=0)
    create_catalogues(maps, 2)

    # populate_cutouts_test((50,50), 5000, (10,15), (5,5))
