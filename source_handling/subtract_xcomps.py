################################################################################
# NAME : subtract_xcomps.py
# DATE STARTED : June 21, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE :
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
# from astropy.FITS_tools.hcongrid import hcongrid #not sure if this is the right module or not, it wasn't clear.
import sys
sys.path.append('../utilities')
from astropy.stats import sigma_clipped_stats
import matplotlib.pyplot as plt
from FITS_tools.hcongrid import hcongrid , hastrom
from get_data import *
from astropy.convolution import convolve
# from scipy import optimize.least_squares
from scipy.stats import linregress
from writefits import *

def subtract_xcomps(maps, simflag=0, verbose=1):
    ncols = len(maps)
    if maps[0]['band'] != 'PSW':
        err = 'first element of map structure is not PSW, aborting.'
        if verbose:
            print(err)
        return None, err


    for i in range(1, ncols):
        if verbose:
            print('On band %s' %(maps[i]['band']))
        # width = sqrt(maps[i]['widtha']**2 - maps[0]['widtha']**2) / maps[0]['pixsize']
        # kern = makeGaussian(15, 15, fwhm =width, center=None)
        # kern = np.array(kern)
        inmap = maps[0]['srcrm']
        # print(inmap.shape, kern.shape)
        print(inmap.shape)
        # inmap = convolve(inmap, kern)
        print(inmap.shape)
        maps[0]['xclean'] = maps[0]['srcrm']
        # whpl = []
        # whnan = []
        # for j in range(inmap.shape[0]):
        #     for k in range(inmap.shape[1]):
        #         maps[0]['mask'][j,k] = 1
        #         inmap[j,k] = 0.0
        # xmap = signal.convolve(inmap, kern)
        xmap = inmap
        print(xmap.shape)
        for j in range(xmap.shape[0]):
            for k in range(xmap.shape[1]):
                pass
                # xmap[WHERE((*maps[0]).mask)] = !VALUES.F_NAN
        xmap_align = hcongrid(xmap, maps[0]['shead'], maps[i]['shead'])

        # mean, median, stddev = sigma_clipped_stats(maps[i]['srcrm'], sigma=10, maxiters=3)


        # whpl = []
        # whnan = []
        # for j in range(maps[i]['srcrm'].shape[0]):
        #     for k in range(maps[i]['srcrm'].shape[1]):
        #         if np.isfinite(maps[i]['srcrm'][j,k]) == True and np.isfinite(xmap[j,k]) == True:
        #             whpl.append([j,k])
        #         else:
        #             whnan.append([j,k])

        #setting up for python version of SVDFIT
        # xmap_align_whpl = []
        # srcrm_whpl = []
        # for value in whpl:
        #     xmap_align_whpl.append(xmap_align[value])
        #     srcrm_whpl.append(maps[i]['srcrm'][value])
        # srcrm_whpl = np.array(srcrm_whpl)
        # xmap_align_whpl = np.array(xmap_align_whpl)

        print(type(maps[i]['srcrm']))
        print(type(xmap_align))
        print(maps[i]['srcrm'])
        for j in range(xmap_align.shape[0]):
            for k in range(xmap_align.shape[1]):
                if np.isnan(xmap_align[j,k]) or np.isnan(maps[i]['srcrm'][j,k]):
                    xmap_align[j,k] = 0
                    maps[i]['srcrm'][j,k] = 0

        print(maps[i]['srcrm'].shape)
        print(xmap_align.shape)
        PSW_array = xmap_align.flatten()
        PMW_array = maps[i]['srcrm'].flatten()

        slope, intercept, r_value, p_value, std_err = linregress(PSW_array, PMW_array)
        y = slope * PSW_array + intercept
        print(r_value)
        # print(np.max(y))
        # print(np.min(y))

        plt.plot(PSW_array, PMW_array, 'x')
        plt.plot(PSW_array, y, c='red')
        plt.show()

        maps[i]['xclean'] = np.empty(maps[i]['srcrm'].shape)

        for j in range(maps[i]['xclean'].shape[0]):
            for k in range(maps[i]['xclean'].shape[1]):
                if maps[i]['srcrm'][j,k] == 0:
                    maps[i]['xclean'][j,k] = np.nan
                else:
                    maps[i]['xclean'][j,k] = maps[i]['srcrm'][j,k] - slope * xmap_align[j,k] + intercept

        # maps[i]['xclean'] = (maps[i]['srcrm'] - slope * xmap_align + intercept)

        datasub = maps[i]['xclean']

        # plt.plot(xmap_align, maps[i]['xclean'], 'x', c='blue')
        # plt.show()
        plt.imshow(maps[i]['xclean'])
        plt.show()

        # x = np.empty(1000)
        # for j in range(x.shape[0]):
        #     x[j] = (j - 500.0) / 1000.0
        # #another call to plot this new thing stil ldon't know if we can just use matplotlib here.
        # for j in range(maps[i]['xclean'].shape[0]):
        #     for k in range(maps[i]['xclean'].shape[1]):
        #         maps[i]['xclean'][j,k] = maps[i]['srcrm'][j,k] - (coeff[1] * xmap_align[j,k] + coeff[0])
        #
        # for indexes in whnan:
        #     maps[i]['mask'][indexes] = 1
        #
        # datasub = maps[i]['xclean']

        if not simflag:
            filename = config.CLUSDATA + 'sz/' + maps[i]['name'] + str(maps[i]['band']) + '_xc.fits'
        else:
            filename = config.CLUSDATA + 'sz/sim/' + maps[i]['name'] + str(maps[i]['band']) + '_xc.fits'

        filename = 'correlated_comp_test_%s.fits' % (maps[i]['band'])
        writefits(filename, data=datasub, header_dict=maps[i]['shead'])
    #
    # for i in range(maps[0]['xclean'].shape[0]):
    #     for j in range(maps[0]['xclean'].shape[1]):
    #         maps[0]['xclean'][i,j] = maps[0]['xclean'][i,j] - mean(maps[0]['xclean'])
    return maps, err

# def makeGaussian(x_size, y_size, fwhm = 3, center=None):
#     """ Make a square gaussian kernel.
#
#     size is the length of a side of the square
#     fwhm is full-width-half-maximum, which
#     can be thought of as an effective radius.
#     """
#
#     sigma = fwhm / 2.355
#     x = np.arange(0, x_size, 1, float)
#     y = np.arange(0, y_size, 1, float)
#     y = y[:,np.newaxis]
#
#     if center is None:
#         x0 = x_size // 2
#         y0 = y_size // 2
#     else:
#         x0 = center[0]
#         y0 = center[1]
#
#     sigma = fwhm / 2.355


if __name__ == '__main__':
    psw = fits.open('../fits_files/xid_9_subtracted_a0370_PSW.fits')
    pmw = fits.open('../fits_files/xid_9_subtracted_a0370_PMW.fits')
    plw = fits.open('../fits_files/xid_9_subtracted_a0370_PLW.fits')
    maps, err = get_data('a0370')
    maps[0]['srcrm'] = psw[0].data
    maps[1]['srcrm'] = pmw[0].data
    maps[2]['srcrm'] = plw[0].data
    print(type(pmw[0].data))
    maps, err = subtract_xcomps(maps)
