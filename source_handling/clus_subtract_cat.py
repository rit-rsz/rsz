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
#   VLB - 04/21/20 - Subtracting off lensing model from PCAT residuals
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
import os, time, math

def clus_subtract_cat(maps, dI, nsim, sgen=None, verbose=1, saveplot=0, superplot=0):
    err = None
    offsets = 3*[0]
    offsets[0] = np.median([x for x in maps[0]['signal'].flatten() if not math.isnan(x)])
    offsets[1] = np.median([x for x in maps[1]['signal'].flatten() if not math.isnan(x)])
    offsets[2] = np.median([x for x in maps[2]['signal'].flatten() if not math.isnan(x)])

    # add 3 mJy offset to make pcat happy for now
    # maps[0]['signal'] = np.reshape([x - offsets[0] + 0.003 for x in maps[0]['signal'].flatten() if not math.isnan(x)],(maps[0]['srcrm'].shape[0],maps[0]['srcrm'].shape[1]))
    # maps[1]['signal'] = np.reshape([x - offsets[1] + 0.003 for x in maps[1]['signal'].flatten() if not math.isnan(x)],(maps[1]['srcrm'].shape[0],maps[1]['srcrm'].shape[1]))
    # maps[2]['signal'] = np.reshape([x - offsets[2] + 0.003 for x in maps[2]['signal'].flatten() if not math.isnan(x)],(maps[2]['srcrm'].shape[0],maps[2]['srcrm'].shape[1]))
    print('map offsets: ', offsets)
    resid_maps = run_pcat(maps,nsim,offsets)

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
            plt.colorbar().set_label('[Jy]')
            plt.title('Clus Subtract Cat : PCAT Residual %s' %(maps[i]['band']))
            plt.savefig(filename)
            plt.clf()

        ''' For making hists for Jack '''
        hda = fits.PrimaryHDU(maps[i]['srcrm'],maps[i]['shead'])
        hda.writeto(config.OUTPUT + 'pcat_residuals/' + maps[i]['name'] + '_resid_' + maps[i]['band'] + '_' + str(nsim) + '.fits',overwrite=True)

        # subtracting lensing template from resids
        data_file = config.HOME + 'Lensing/lense_template_' + maps[i]['band'] + '.fits'
        lense_model = fits.getdata(data_file)
        maps[i]['srcrm'] = maps[i]['srcrm'] - lense_model

    return maps, err

def run_pcat(maps,nsim,offsets):
    ob = lion(band0=0, band1=1, band2=2, mean_offsets = offsets, isim=nsim, map_object=maps, auto_resize=True, make_post_plots=True, openblas=False, cblas=False, nsamp=500, residual_samples=100, visual=False)
    resid_maps = ob.main()
    return resid_maps
