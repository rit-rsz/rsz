################################################################################
# NAME : fitsz.py
# DATE STARTED : June 24, 2019
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
import numpy as np
import matplotlib.pyplot as plt
import sys
from math import *
sys.path.append('../utilities')
from get_clusparams import *
from compute_rings import *
from get_data import *
from scipy.stats import linregress


def fitsz(radave, params, beam=None, maxlim=3600, minlim=0, noweight=1, superplot=1, verbose=1):

    # init params
    ncols = len(radave)
    increment = np.zeros(ncols)
    offsets = np.zeros(ncols)
    fit = np.zeros((2, ncols))

    for i in range(ncols):
        print(params['rc'])
        print(radave[i]['midbin'])
        print(radave[i]['midbin'].shape)
        rc = np.tile(params['rc'], radave[i]['midbin'].shape)
        beta = np.tile(params['beta'], radave[i]['midbin'].shape)
        if  not beam:
            roft = (1. + (radave[i]['midbin'] / rc)**2)**((1. - 3. * beta / 2.))

        #     for j in range(radave[i]['midbin'].shape[0]):
        #         for k in range(radave[i]['midbin'].shape[1]):
        #             roft = (1.0 + (radave[i]['midbin'] / rc)**2)**((1.0 - 3.0 * beta) / 2.0)
        else:
            print('this part is not written yet')
            #this stuff is commented out because I don't understand it also it's probably wrong.
            # xrad = []
            # for j in range(np.amax.(radave[i]['midbin']), 0, -1):
            #     xrad.append(j)
            # xrad = np.array(xrad)
            # xrc = np.tile(params['rc'],2*np.amax(radave[i]['midbin']))
            # xbeta = np.tile(params['beta'], 2*np.amax(radave[i]['midbin']))
            # roftprime = (1. + (xrad / xrc)**2)**((1. - 3. * xbeta) / 2.)
            # pf = -(2.18 * log(2.) / beam)
            # psf = np.empty(len(xrad))
            # for j in range(len(radave[i]['midbin'])):
            #     roftprime[j] = (1.0 + (xrad[j] / xrc[j])**2)**((1.0 - 3.0 * xbeta[j]) / 2.0)
            #     pf = -1*(2.18 * math.log(2.0) / beam[i])**2
            #     psf[j] = math.exp(pf * xrad[j]**2)
            #     psf[j] = psf[j] / np.sum(psf)
            #     #roftprimep = CONVOL(roftprime, psf, EDGE_TRUNCATE)
            #     roftprimep[j] = roftprimep[j] / np.aload_headermax(abs(roftprimep))
            #     roft = INTERPL(roftprimep, xrad, radwave[i]['midbin'])

        if superplot:
            plt.plot(radave[i]['midbin'], radave[i]['fluxbin'])
            plt.show()
            #code to plot the stuffs this could probably be done in matplotlib...

        ind = []
        for j in range(len(roft)):
            if np.isfinite(radave[i]['fluxbin'][j]) and radave[i]['midbin'][j] <= maxlim:
                ind.append(j)
        ind = np.asarray(ind)

        slope, intercept, r_value, p_value, std_err = linregress(roft[ind], radave[i]['fluxbin'][ind])

        A = slope
        B = intercept

        print(A, B)

        if superplot:
            y = slope * roft + intercept
            plt.plot(roft, radave[i]['fluxbin'], c='blue')
            plt.plot(roft, y, c='red')
            plt.show()
            #plot the new stuffs, again could probably be done in matplotlib, but i need to do real testing to find out.

        increment[i] = slope
        offsets[i] = intercept
        fit[:, i] = np.array([slope, intercept])

    return fit, offsets, increment

if __name__ == '__main__':
    params, err = get_clus_params('a0370')
    maps, err = get_data('a0370')
    psw = fits.open('../fits_files/xid_9_subtracted_a0370_PSW.fits')
    pmw = fits.open('../fits_files/xid_9_subtracted_a0370_PMW.fits')
    plw = fits.open('../fits_files/xid_9_subtracted_a0370_PLW.fits')
    maps[0]['srcrm'] = psw[0].data
    maps[1]['srcrm'] = pmw[0].data
    maps[2]['srcrm'] = plw[0].data
    radave = compute_rings(maps,params,30.0)
    print(radave[0]['midbin'])
    fit = fitsz(radave, params)
