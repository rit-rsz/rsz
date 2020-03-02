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
from math import *
import sys
sys.path.append('../utilities')
import matplotlib.pyplot as plt
from FITS_tools.hcongrid import hcongrid , hastrom
from clus_get_data import *
from astropy.convolution import convolve_fft
from scipy.stats import linregress
from save_fits import writefits
from gaussian import makeGaussian, padGaussian
from clus_make_noise_mask import clus_make_noise_mask
import config

def clus_subtract_xcomps(maps, sgen=None, verbose=1, superplot=1, nsim=0, saveplot=1):

    err =  None
    ncols = len(maps)
    #check to see if the maps are in the correct order.
    if maps[0]['band'] != 'PSW':
        err = 'first element of map structure is not PSW, aborting.'
        if verbose:
            print(err)
        return None, err

    #make a noise mask and populate source removed map with nans where mask is 1 so they don't effect the fit.
    for i in range(len(maps)):
        mask = clus_make_noise_mask(maps, i)
        maps[i]['mask'] = maps[i]['mask'] + mask
        for j in range(maps[i]['signal'].shape[0]):
            for k in range(maps[i]['signal'].shape[1]):
                if maps[i]['mask'][j,k] == 1:
                    maps[i]['srcrm'][j,k] = np.nan

        if saveplot :
            if sgen != None :
                filename = config.OUTPUT + 'pcat_residuals/' + maps[i]['name'] + '_mask_resid_' + maps[i]['band'] + '_' + str(nsim) + '.png'
            else :
                filename = config.OUTPUT + 'pcat_residuals/' + maps[i]['name'] + '_mask_resid_' + maps[i]['band'] + '_real' + '.png'
            plt.imshow(maps[i]['srcrm'])
            plt.colorbar().set_label('[Jy]')
            plt.title('Clus Subtract Xcomps : Masked PCAT Residual %s' %(maps[i]['band']))
            plt.savefig(filename)
            plt.clf()

    for i in range(1, ncols):
        if verbose:
            print('On band %s' %(maps[i]['band']))

        #create a new image for the PSW band that is the same shape and beamsize as the reference images
        width = sqrt(maps[i]['widtha']**2 - maps[0]['widtha']**2) / maps[0]['pixsize']
        kern = makeGaussian(15, 15, fwhm =width, center=None)
        kern = kern / np.sum(kern)
        kern = padGaussian(maps[0]['srcrm'],kern)
        inmap = convolve_fft(maps[0]['srcrm'], kern)

        inmap = maps[0]['srcrm']
        maps[0]['xclean'] = maps[0]['srcrm']
        xmap = inmap
        xmap_align = hcongrid(xmap, maps[0]['shead'], maps[i]['shead'])

        # plt.imshow(xmap_align)
        # plt.colorbar()
        # plt.clim(-0.03,0.04)
        # plt.savefig('xmap_align_%s.png' %(maps[i]['band']))
        # plt.clf()

        for j in range(xmap_align.shape[0]):
            for k in range(xmap_align.shape[1]):
                if np.isnan(xmap_align[j,k]) or np.isnan(maps[i]['srcrm'][j,k]):
                    xmap_align[j,k] = 0
                    maps[i]['srcrm'][j,k] = 0

        #now that we have our new PSW image flatten both images
        PSW_array = xmap_align.flatten()
        ref_array = maps[i]['srcrm'].flatten()

        #find the linear fit for PSW vs ref
        slope, intercept, r_value, p_value, std_err = linregress(PSW_array, ref_array)
        y = slope * PSW_array + intercept


        if superplot or saveplot:
            plt.plot(PSW_array, ref_array, 'x',label='Raw Intensities')
            plt.plot(PSW_array, y, c='red', label='Linear Regression Fit')
            plt.legend()
            plt.title('Clus Subtract Xcomps : PSW vs. %s' %(maps[i]['band']))
            plt.xlabel('PSW [Jy]')
            plt.ylabel('%s [Jy]' %(maps[i]['band']))

            if superplot:
                plt.show()
            elif saveplot:
                if sgen != None:
                    filename = config.OUTPUT + 'corr_comps/' + maps[i]['name'] + '_xcomps_' + maps[i]['band'] + '_' + str(nsim) + '.png'
                else:
                    filename = config.OUTPUT + 'corr_comps/' + maps[i]['name'] + '_xcomps_' + maps[i]['band'] + '_real.png'

                plt.savefig(filename)
                plt.clf()

        #subtract the correlated components from the image
        maps[i]['xclean'] = np.empty(maps[i]['srcrm'].shape)
        for j in range(maps[i]['xclean'].shape[0]):
            for k in range(maps[i]['xclean'].shape[1]):
                if maps[i]['srcrm'][j,k] == 0:
                    maps[i]['xclean'][j,k] = np.nan
                else:
                    maps[i]['xclean'][j,k] = maps[i]['srcrm'][j,k] - slope * xmap_align[j,k] + intercept

        plt.imshow(maps[i]['xclean'])
        plt.colorbar().set_label('[Jy]')
        plt.title('Clus Subtract Xcomp : Correlated Components Removed')

        if superplot:
            plt.show()
        elif saveplot:
            if sgen != None:
                filename = config.OUTPUT + 'corr_comps/' + maps[i]['name'] + '_xclean_' + maps[i]['band'] + '_' + str(nsim) + '.png'
            else:
                filename = config.OUTPUT + 'corr_comps/' + maps[i]['name'] + '_xclean_' + maps[i]['band'] + '_real.png'

            plt.savefig(filename)
            plt.clf()

    #subtract the mean of the new map from itself. Why do we do this? we need to figure out what the purpose of this step is.
    # for i in range(maps[0]['xclean'].shape[0]):
    #     for j in range(maps[0]['xclean'].shape[1]):
    #         maps[0]['xclean'][i,j] = maps[0]['xclean'][i,j] - np.mean(maps[0]['xclean'])

    return maps, err
