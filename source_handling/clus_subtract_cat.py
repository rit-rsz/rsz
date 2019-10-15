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
import matplotlib.pyplot as plt
# from clus_get_data import *
from astropy.io import fits
import config

def clus_subtract_cat(maps, cat, verbose=1):
    err = False

    # make sure both input maps exist
    if maps == None or cat == None :
        return None, 'clus or pcat_spire map structure input is absent'

    for i in range(len(cat)): # only loop through how many residuals we have
        if verbose:
            print('Subtracting for %s' %(maps[i]['band']))

        # make the difference image
        datasub = np.empty(maps[i]['astr']['NAXIS'])
        # make full data map object
        datafull = maps[i]['signal'].data

        # find parts of the signal map & mask which exist and dont
        # perform subtraction for good data, make nans for bad data
        for j in range(maps[i]['signal'].data.shape[0]):
            for k in range(maps[i]['signal'].data.shape[1]):
                if np.isfinite(maps[i]['signal'].data[j,k]) == False or maps[i]['mask'][j,k] == 1 :
                    datasub[j,k] = np.nan
                else:
                    datasub[j,k] = maps[i]['signal'].data[j,k] - cat[i][j,k]

        # write a fits for debugging purposes
        # if os.path.exists(config.FITSOUT + 'xid_subtracted_%s_%s.fits' % (maps[i]['name'], maps[i]['band'])) == True :
        #     os.remove(config.FITSOUT + 'xid_subtracted_%s_%s.fits' % (maps[i]['name'], maps[i]['band']))
        #     hdu = fits.PrimaryHDU(datasub)
        #     hdul = fits.HDUList([hdu])
        #     hdul.writeto(config.FITSOUT + 'xid_subtracted_%s_%s.fits' % (maps[i]['name'], maps[i]['band']))

        # plotting for debugging purposes
        fig1, ax1 = plt.subplots()
        fig1, ax2 = plt.subplots()
        cp1 = ax1.contour(datafull)
        ax1.set_title('clus_subtract_cat: Signal map for %s' %(maps[i]['band']))
        cp2 = ax2.contour(datasub)
        ax2.set_title('clus_subtract_cat: Catalog subtracted map for %s' %(maps[i]['band']))
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
    resid_maps = np.load('/home/butler/rsz/pcat_resid.npy')
    subtract_cat(maps,resid_maps)
