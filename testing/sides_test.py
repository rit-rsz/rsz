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

def map_selector(m_size, s_size, n_samp):
    """
    m_size : catalog size
    s_size : cut out size
    n_samp : number of cuttouts
    """
    #create a new map that allows us to select center pixels for potential SPIRE
    #maps that don't overextend the edges of the SIDES map.

    #calculate needed floating point values to generate centers within the correct
    #area
    n_shape = (m_size[0] - floor(s_size[0]), m_size[1] - 2*floor(s_size[1] / 2))
    #generate random center pixels with x, y coordinates
    ### testing code
    c_pix = []
    for i in range(n_samp):
        cx = np.random.random() * n_shape[0] + floor(s_size[0] / 2)
        cy = np.random.random() * n_shape[1] + floor(s_size[1] / 2)
        #check the newly generated coordinates with a list of previously generated coordinates
        c = check_c(c_pix, [cx, cy], n_shape, s_size)
        #plotting squares for testing purposes
        plot_squares([c[0], c[1]], s_size)
        #append the center pixel to our list
        c_pix.append(c)
    c_pix = np.asarray(c_pix)

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


    plt.plot(ix, iy)
    plt.scatter(c_pix[:,0], c_pix[:,1])

    plt.plot(gx, gy, c='red')
    plt.show()

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


def calc_ds(c_list, c_pix):
    #calculate manhattan distance and return a list of distances
    d_list = []
    for c in c_list:
        d = abs(c[0] - c_pix[0]) + abs(c[1] - c_pix[1])
        if d < 25:
            d_list.append(d)
    return d_list

def check_c(c_list, c_pix, n_shape, s_size):
    iter = 0
    if len(c_list) > 0:
        d_list = calc_ds(c_list, c_pix)
        while len(d_list) > 0:
            if iter > 100:
                print('Error: exceeded iteration limit.')
                break #break the while loop if it goes over 100 iterations
            iter += 1
            cx = np.random.random() * n_shape[0] + floor(s_size[0] / 2)
            cy = np.random.random() * n_shape[1] + floor(s_size[1] / 2)
            d_list = calc_ds(c_list, [cx, cy])
        return c_pix
    else:
        return c_pix

if __name__ == '__main__':
    map_selector((1000,1000), (100,100), 100)
