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
from astropy.io import fits
sys.path.append('/utilities')
import config
sys.path.append(config.HOME + 'multiband_pcat/multiband_pcat')
from pcat_spire import lion
from save_fits import writefits
import os, time

def clus_subtract_cat(maps, dI, nsim, sgen=None, verbose=1, saveplot=0, superplot=0):
    err = None

    resid_maps = run_pcat(maps)

    for i in range(len(resid_maps)): # only loop through how many residuals we have
        if verbose==1:
            print('Setting Subtracted Maps for %s' %(maps[i]['band']))

        # make the difference image
        datasub = resid_maps[i]
        # make full data map object
        datafull = maps[i]['signal']

        if superplot:
            # plotting for debugging purposes
            fig1, ax1 = plt.subplots()
            fig1, ax2 = plt.subplots()
            cp1 = ax1.contour(datafull)
            ax1.set_title('clus_subtract_cat: Signal map for %s' %(maps[i]['band']))
            cp2 = ax2.contour(datasub)
            ax2.set_title('clus_subtract_cat: Catalog subtracted map for %s' %(maps[i]['band']))
            plt.show()

            plt.imshow(datasub)
            plt.colorbar()
            plt.title('clus_subtract_cat: Catalog subtracted map for %s' %(maps[i]['band']))
            plt.show()

        # save new subtracted signal to map
        # reshape pcat square map to original SPIRE size
        maps[i]['srcrm'] = datasub[0:maps[i]['signal'].shape[0],0:maps[i]['signal'].shape[1]]

        if saveplot:
            if sgen != None:
                filename = config.OUTPUT + 'pcat_residuals/' + maps[i]['name'] + '_resid_' + maps[i]['band'] + '_' + str(nsim) + '.png'
            else:
                filename = config.OUTPUT + 'pcat_residuals/' + maps[i]['name'] + '_resid_' + maps[i]['band'] + '_' + 'real.png'

            plt.imshow(maps[i]['srcrm'])
            plt.colorbar()
            plt.savefig(filename)
            plt.clf()


    return maps, err

def run_pcat(maps):
    ob = lion(band0=0, band1=1, band2=2, map_object=maps, auto_resize=True, make_post_plots=True, openblas=False, cblas=False, nsamp=1, residual_samples=1, visual=False)
    resid_maps = ob.main()
    return resid_maps
