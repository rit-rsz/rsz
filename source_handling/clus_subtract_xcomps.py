################################################################################
# NAME : subtract_xcomps.py
# DATE STARTED : June 21, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE : The purpose of this routine is to subtract the correlated components
# of a reference image with that of the PSW image
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
from math import *
from scipy import signal
from math import *
import sys
sys.path.append('../utilities')
from astropy.stats import sigma_clipped_stats
import matplotlib.pyplot as plt
from FITS_tools.hcongrid import hcongrid , hastrom
from clus_get_data import *
from astropy.convolution import convolve_fft
# from scipy import optimize.least_squares
from scipy.stats import linregress
from writefits import *
from gaussian import makeGaussian, padGaussian

def clus_subtract_xcomps(maps, simflag=0, verbose=1, superplot=1):

    err =  None
    ncols = len(maps)
    #check to see if the maps are in the correct order.
    if maps[0]['band'] != 'PSW':
        err = 'first element of map structure is not PSW, aborting.'
        if verbose:
            print(err)
        return None, err

    mask = clus_make_noise_mask(maps, 0)
    maps[0]['mask'] = maps[0]['mask'] + mask

    for i in range(1, ncols):
        if verbose:
            print('On band %s' %(maps[i]['band']))

        #create a new image for the PSW band that is the same shape and beamsize as the reference images
        width = sqrt(maps[i]['widtha']**2 + maps[0]['widtha']**2) / maps[0]['pixsize']
        kern = makeGaussian(15, 15, fwhm =width, center=None)
        inmap = maps[0]['srcrm']
        kern = padGaussian(inmap,kern)
        inmap = convolve_fft(inmap, kern)
        maps[0]['xclean'] = maps[0]['srcrm']
        xmap = inmap
        xmap_align = hcongrid(xmap, maps[0]['shead'], maps[i]['shead'])

        mask = clus_make_noise_mask(maps, i)

        for j in range(xmap_align.shape[0]):
            for k in range(xmap_align.shape[1]):
                maps[m]['mask'][ipix,jpix] = maps[m]['mask'][ipix,jpix] + mask[ipix,jpix]
                if np.isnan(xmap_align[j,k]) or np.isnan(maps[i]['srcrm'][j,k]) or maps[m]['mask'][ipix,jpix] == 1:
                    xmap_align[j,k] = 0
                    maps[i]['srcrm'][j,k] = 0

        #now that we have our new PSW image flatten both images
        PSW_array = xmap_align.flatten()
        ref_array = maps[i]['srcrm'].flatten()

        #find the linear fit for PSW vs ref
        slope, intercept, r_value, p_value, std_err = linregress(PSW_array, ref_array)
        y = slope * PSW_array + intercept


        if superplot:
            plt.plot(PSW_array, ref_array, 'x')
            plt.plot(PSW_array, y, c='red')
            plt.title('clus_subtract_xcomps: PSW vs. %s' %(maps[i]['band']))
            plt.xlabel('PSW')
            plt.ylabel('%s' %(maps[i]['band']))
            plt.show()

        #subtract the correlated components from the image
        maps[i]['xclean'] = np.empty(maps[i]['srcrm'].shape)
        for j in range(maps[i]['xclean'].shape[0]):
            for k in range(maps[i]['xclean'].shape[1]):
                if maps[i]['srcrm'][j,k] == 0:
                    maps[i]['xclean'][j,k] = np.nan
                else:
                    maps[i]['xclean'][j,k] = maps[i]['srcrm'][j,k] - slope * xmap_align[j,k] + intercept

        datasub = maps[i]['xclean']

        if superplot:
            plt.imshow(maps[i]['xclean'])
            plt.show()

        if not simflag:
            filename = config.CLUSDATA + 'sz/' + maps[i]['name'] + str(maps[i]['band']) + '_xc.fits'
        else:
            filename = config.CLUSDATA + 'sz/sim/' + maps[i]['name'] + str(maps[i]['band']) + '_xc.fits'

        filename = 'correlated_comp_test_%s.fits' % (maps[i]['band'])
        if os.path.isfile(filename):
            os.remove(filename)
        writefits(filename, data=datasub, header_dict=maps[i]['shead'])

    #subtract the mean of the new map from itself.
    for i in range(maps[0]['xclean'].shape[0]):
        for j in range(maps[0]['xclean'].shape[1]):
            maps[0]['xclean'][i,j] = maps[0]['xclean'][i,j] - np.mean(maps[0]['xclean'])

    plt.imshow(maps[0]['xclean'])
    plt.show()

    return maps, err



if __name__ == '__main__':
    from clus_get_data import clus_get_data
    maps, err = clus_get_data('a0370')

    # psw = fits.open('../fits_files/xid_9_subtracted_a0370_PSW.fits')
    # pmw = fits.open('../fits_files/xid_9_subtracted_a0370_PMW.fits')
    # plw = fits.open('../fits_files/xid_9_subtracted_a0370_PLW.fits')
    # maps[0]['srcrm'] = psw[0].data
    # maps[1]['srcrm'] = pmw[0].data
    # maps[2]['srcrm'] = plw[0].data

    import sys
    sys.path.append('../utilities')
    sys.path.append('../source_handling')
    import numpy as np
    map = np.load('/home/butler/rsz/pcat_resid.npy',allow_pickle=True)
    map1 = np.load('/home/butler/rsz/pcat_resid1.npy',allow_pickle=True)
    map2 = np.load('/home/butler/rsz/pcat_resid2.npy',allow_pickle=True)
    maps[0]['srcrm'] = map[0][0:maps[0]['signal'].shape[0],0:maps[0]['signal'].shape[1]]
    maps[1]['srcrm'] = map1[0][0:maps[1]['signal'].shape[0],0:maps[1]['signal'].shape[1]]
    maps[2]['srcrm'] = map2[0][0:maps[2]['signal'].shape[0],0:maps[2]['signal'].shape[1]]
    from clus_subtract_xcomps import clus_subtract_xcomps
    clus_subtract_xcomps(maps)
