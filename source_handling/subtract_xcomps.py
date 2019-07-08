################################################################################
# NAME : subtract_xcomps.py
# DATE STARTED : June 21, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE : idk lol
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
from meanclip import meanclip

def subtract_xcomps(maps, simflag=0, verbose=1):
    ncols = len(maps)

    if maps[0]['band'] != 'PSW':
        err = 'first element of map structure is not PSW, aborting.'
        if verbose:
            print(err)
        return None, err

    for i in range(ncols):
        if verbose:
            print('On band %s' %(maps[i]['band']))
        width = sqrt(maps[i]['widtha']**2 - maps[0]['widtha']**2) / maps[0]['pixsize']
        stdev = width / (2*sqrt(2 * log(2.)))
        kern = Gaussian2DKernel(stdev, x_size=15, y_size=15)
        kern = np.array(kern)
        inmap = maps[0]['srcrm']
        maps[0]['xclean'] = maps[0]['srcrm']
        whpl = []
        whnan = []
        for j in range(inmap.shape[0]):
            for k in range(inmap.shape[1]):
                if np.isfinite(inmap[j,k]) == True:
                    whpl.append([j,k])
                else:
                    whnan.append([j,k])
        for indexes in whnan:
            maps[0]['mask'][indexes] = 1
            inmap[indexes] = 0.0
        xmap = signal.convolve2d(inmap, kern)
        indexes = []
        for j in range(xmap.shape[0]):
            for k in range(xmap.shape[1]):
                if maps[0]['mask'][j,k] != 0:
                    indexes.append([j,k])
        for value in indexes:
            xmap[value] = np.nan #not sure about this because it doesn't specify a condition for where in IDL so i don't know if my if statement is correct.
        xmap_align = hcongrid(xmap, maps[0]['shead'], maps[i]['shead']) #this may not work the package was FITS_tools.hcongrid.hcongrid but I'm not sure if its part of astropy or not.

        mean, sigma = meanclip(maps[i]['srcrm'], clipsig=10, maxiter=3, verbose=verbose)
        whpl = []
        whnan = []
        for j in range(maps[i]['srcrm'].shape[0]):
            for k in range(maps[i]['srcrm'].shape[1]):
                if np.isfinite(maps[icol]['srcrm'][j,k]) == True and np.isfinite(xmap[j,k]) == True:
                    whpl.append([j,k])
                else:
                    whnan.append([j,k])
        #then it calls plot, so i don't know if we can just matplotlib here but ok i guess

        #setting up for python version of SVDFIT
        xmap_align_whpl = []
        srcrm_whpl = []
        for value in whpl:
            xmap_align_whpl.append(xmap_align[value])
            srcrm_whpl.append(maps[i]['srcrm'][value])
        srcrm_whpl = np.array(srcrm_whpl)
        xmap_align_whpl = np.array(xmap_align_whpl)

        #call to SVDFIT which is a least squares fit.


        x = np.empty(1000)
        for j in range(x.shape[0]):
            x[j] = (j - 500.0) / 1000.0
        #another call to plot this new thing stil ldon't know if we can just use matplotlib here.
        for j in range(maps[i]['xclean'].shape[0]):
            for k in range(maps[i]['xclean'].shape[1]):
                maps[i]['xclean'][j,k] = maps[i]['srcrm'][j,k] - (coeff[1] * xmap_align[j,k] + coeff[0])

        for indexes in whnan:
            maps[i]['mask'][indexes] = 1

        datasub = maps[i]['xclean']

        if not simflag:
            filename = config.CLUSDATA + 'sz/' + maps[i]['name'] + str(maps[i]['band']) + '_xc.fits'
        else:
            filename = config.CLUSDATA + 'sz/sim/' + maps[i]['name'] + str(maps[i]['band']) + '_xc.fits'
        writefits(filename, data=datasub, header_dict=maps[i]['shead'])
        #another call to contour need to fix this.

    for i in range(maps[0]['xclean'].shape[0]):
        for j in range(maps[0]['xclean'].shape[1]):
            maps[0]['xclean'][i,j] = maps[0]['xclean'][i,j] - mean(maps[0]['xclean']) # this is questionable for a 2d array.

    return maps, err
