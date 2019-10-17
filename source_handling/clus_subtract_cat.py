################################################################################
# NAME : subtract_cat.py
# DATE STARTED : June 21, 2019
# AUTHORS : Benjamin Vaughan & Victoria Butler
# PURPOSE : the purpose of this script is to subtract the data from XID from
#           the original signal
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
#   VLB - 8/13/19 - trimmed down the code to get rid of unneccessary for loops,
#                   added contour plots to visualize output, robust deletion and
#                   creation of new fits files for debugging
#   VLB - 10/15/19 - changing to work with PCAT_SPIRE residual maps instead of XID+
################################################################################
import numpy as np
import sys
import matplotlib.pyplot as plt
# from clus_get_data import *
from astropy.io import fits
sys.path.append('../utilities')
import config
sys.path.append(config.HOME + 'multiband_pcat')
from pcat_spire import lion

def clus_subtract_cat(maps, verbose=1):
    err = None

    # ob = lion(raw_counts=True, auto_resize=True, visual=True)
    # default is to have all of the 3 maps returned
    ob = lion(map_object=maps, auto_resize=True, make_post_plots=False, nsamp=100, residual_samples=100, visual=False)
    resid_maps = ob.main()

    plt.imsave('/home/butler/rsz/pcat_resid.png',resid_maps[0],format='png')
    np.save('/home/butler/rsz/pcat_resid.npy',resid_maps)
    # resid_maps = np.load('/home/butler/rsz/pcat_resid.npy',allow_pickle=True)

    # make sure both input maps exist
    ''' turned off for testing purposes '''
    # if maps.any() == None or resid_maps.any() == None :
    #     return None, 'clus or pcat_spire map structure input is absent'

    for i in range(len(resid_maps)): # only loop through how many residuals we have
        if verbose==1:
            print('Setting Subtracted Maps for %s' %(maps[i]['band']))

        # make the difference image
        datasub = resid_maps[i]
        # make full data map object
        datafull = maps[i]['signal']

        # plotting for debugging purposes
        # fig1, ax1 = plt.subplots()
        # fig1, ax2 = plt.subplots()
        # cp1 = ax1.contour(datafull)
        # ax1.set_title('clus_subtract_cat: Signal map for %s' %(maps[i]['band']))
        # cp2 = ax2.contour(datasub)
        # ax2.set_title('clus_subtract_cat: Catalog subtracted map for %s' %(maps[i]['band']))
        # plt.show()

        plt.imshow(datasub)
        plt.colorbar()
        plt.title('clus_subtract_cat: Catalog subtracted map for %s' %(maps[i]['band']))
        plt.show()
        # save new subtracted signal to map
        maps[i]['scrrm'] = datasub

    return maps, err

if __name__ == '__main__' :
    import sys
    sys.path.append('../utilities')
    sys.path.append('../source_handling')
    from clus_get_data import clus_get_data
    maps, err = clus_get_data(clusname='a0370')
    resid_maps = np.load('/home/butler/rsz/pcat_resid.npy',allow_pickle=True)
    clus_subtract_cat(maps,resid_maps)
