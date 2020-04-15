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
    min_pix = np.min(pset) #find max and min pixel values for limits of histogram
    max_pix = np.max(pset)
    N = round((max_pix - min_pix) / std)
    bins = np.linspace(min_pix, max_pix, N) #create bins based on the uncertainty in pixel vals
    counts, _ = np.histogram(pset, bins)
    print(len(counts), len(bins[:-1]))
    plt.scatter(bins[:-1], counts)
    plt.show()
    c, b, _ = plt.hist(pset, N, type='step') #plotting for debugging
    plt.show()


if __name__ == '__main__':
    sim_one_pixel(3, 1) #simulate with background of 3 and std of 1
