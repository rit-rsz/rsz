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
    #create a new map that allows us to select center pixels for potential SPIRE
    #maps that don't overextend the edges of the SIDES map.
    n_shape = (m_size[0] - floor(s_size[0] / 2), m_size[1] - floor(s_size[1] / 2))
    c_pix = []
    for i in range(n_samp):
        cx = np.random.random() * n_shape[0]
        cy = np.random.random() * n_shape[1]
        plot_squares([cx, cy], s_size)
        c_pix.append([cx, cy])

    #plot the greater map
    gc1 = [0,0]
    gc2 = [0,m_size[1]]
    gc3 = [m_size[0], m_size[1]]
    gc4 = [m_size[0], 0]
    ggrid = [gc1, gc2, gc3, gc4, gc1]
    plt.plot(ggrid)
    plt.show()
    #show that the centers are within an area such that edges of the map don't
    #overextend the SIDES catalog
    #find the inner corners
    iup_x = 0 + s_size[0] / 2.
    iup_y = 0 + s_size[1] / 2.
    idown_x = m_size[0] - s_size[0] / 2.
    idown_y = m_size[1] - s_size[1] / 2.

    ic1 = [idown_x, idown_y]
    ic2 = [idown_x, iup_y]
    ic3 = [iup_x, iup_y]
    ic4 = [iup_x, idown_y]
    igrid = [ic1, ic2, ic3, ic4, ic1]
    plt.plot(igrid)
    plt.scatter(c_pix)
    plt.show()





def plot_squares(c_pix, s_size):
    #break up inputs into individual values
    cx = c_pix[0]
    cy = c_pix[1]
    sx_2 = s_size[0] / 2.
    sy_2 = s_size[1] / 2.
    #find the corner pieces of the square
    up_x = cx + sx_2
    down_x = cx - sx_2
    up_y = cy + sy_2
    down_y = cy - sy_2
    c1 = [up_x, down_y]
    c2 = [up_x, up_y]
    c3 = [down_x, up_y]
    c4 = [down_x, down_y]
    #create a list of points to draw the grid
    grid = [c1, c2, c3, c4, c1]
    #plot the grid
    plt.plot(grid)

if __name__ == '__main__':
    map_selector((1000,1000), (100,100), 5)
