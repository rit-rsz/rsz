############################################################################
# NAME : clus_make_lensed_bkgrd.py
# DATE : April 14, 2020
# AUTHOR : Ben Vaughan
# PURPOSE : The purpose of this routine is to create a mean reaization of the
# lensed background after lensing has been applied to the map.
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
# OUTPUTS :
# REVISION HISTORY :

############################################################################
import numpy as np
import matplotlib.pyplot as plt
from math import *

def create_bin_edges(pset, std):
    #the purpose of this function is to create the binning for our pixels
    '''
    pset = the set of pixels at a given coordinate in each image
    std  = the std of the average of the pixels
    '''
    min_pix = np.min(pset) #find max and min pixel values for limits of histogram
    max_pix = np.max(pset)
    N = ceil((max_pix - min_pix) / std)#omething to note is that the number of bins has to be an integer
    #this causes an associated error in the binning that has to be addressed.
    stop = min_pix + std * N #this is a potential solution to address the stepping problem. (so far does not work)
    bin_cent, stepcheck = np.linspace(min_pix, stop, N+1, retstep=True) #create bins based on the uncertainty in pixel vals
    #something to note is that the bin edges correspond to the left side in general + the last right edge
    # for example binning integers 1 - 4we would have [1,2) for 1-2 and [2,3) for 2 -3 and [3,4] for 3 - 4.
    left_edges = np.add(bin_cent, -0.5 * std) #the edges left justified edges of our  bins
    right_edge = np.array(left_edges[-1] + std) #the last edge which should be right justified
    bin_edges = np.append(left_edges, right_edge    ) #add the last right justified bin to get plotting to work
    return bin_edges, bin_cent

def sim_one_pixel(base, std):
    #the purpose of this function is to simulate the creation of a histogram from
    '''
    base = base level of background
    std  = square root of the variance of the background
    '''
    #100 different maps
    noise = np.random.normal(0, std, 100) #create gaussian noise centered around 0 with a known sigma
    pset = np.add(noise, base) #set of pixels with some base background with gaussian noise
    #create a histogram
    bins, cent = create_bin_edges(pset, std)
    counts, bincheck1 = np.histogram(pset, bins)
    #plots for debugging
    countcheck, bincheck2, _ = plt.hist(pset, bins, histtype='step')
    plt.scatter(cent, counts)
    plt.show()


if __name__ == '__main__':
    sim_one_pixel(3, 1) #simulate with background of 3 and std of 1
