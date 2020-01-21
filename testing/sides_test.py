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
    cx = np.random.random() * n_shape[0]
    cy = np.random.random() * n_shape[1]

def plot_squares(c_pix):
               

if __name__ == '__main__':
    map_selector((1000,1000), (100,100), 5)
