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
from save_fits import writefits
import os

def clus_subtract_cat(maps, verbose=1, nsim=0, saveplot=1):
    err = None

    # default is to have all of the 3 maps returned
    ob = lion(band0=0, band1=1, band2=2, map_object=maps, auto_resize=True, make_post_plots=False, nsamp=500, residual_samples=100, visual=False)
    resid_maps = ob.main()

    # plt.imsave('/home/butler/rsz/pcat_resid_0.png',resid_maps[0],format='png')
    # plt.imsave('/home/butler/rsz/pcat_resid_1.png',resid_maps[1],format='png')
    # plt.imsave('/home/butler/rsz/pcat_resid_2.png',resid_maps[2],format='png')
    # np.save('/home/butler/rsz/pcat_resid.npy',resid_maps)

    # map1 = np.load('/home/butler/rsz/pcat_resid.npy',allow_pickle=True)
    # map2 = np.load('/home/butler/rsz/pcat_resid1.npy',allow_pickle=True)
    # map3 = np.load('/home/butler/rsz/pcat_resid2.npy',allow_pickle=True)

    # resid_maps = [map1[0],map2[0],map3[0]]

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

        # plt.imshow(datasub)
        # plt.colorbar()
        # plt.title('clus_subtract_cat: Catalog subtracted map for %s' %(maps[i]['band']))
        # plt.show()
        # save new subtracted signal to map
        # reshape pcat square map to original SPIRE size
        maps[i]['srcrm'] = datasub[0:maps[i]['signal'].shape[0],0:maps[i]['signal'].shape[1]]
        if saveplot:
            if nsim != 0:
                filename = config.HOME + 'outputs/pcat_residuals/' + maps[i]['name'] + '_' + maps[i]['band'] + '_' + str(nsim) + '.fits'
            else:
                filename = config.HOME + 'outputs/pcat_residuals/' + maps[i]['name'] + '_' + maps[i]['band'] + '_' + '.fits'
            if os.path.isfile(filename):
                os.remove(filename)
            writefits(filename, data=maps[i]['srcrm'])


    return maps, err



if __name__ == '__main__' :
    import sys
    sys.path.append('../utilities')
    sys.path.append('../source_handling')
    from clus_get_data import clus_get_data
    maps, err = clus_get_data(clusname='a0370')
    resid_maps = np.load('/home/butler/rsz/pcat_resid.npy',allow_pickle=True)
    clus_subtract_cat(maps,resid_maps)
