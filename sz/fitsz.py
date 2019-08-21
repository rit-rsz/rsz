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
from astropy.convolution import convolve_fft
from scipy.interpolate import splrep
from gaussian import *
from get_spire_beam_fwhm import *

def fitsz(radave, params, beam=None, maxlim=3600, minlim=0, noweight=1, superplot=1, verbose=1):

    # init params
    beam = [get_spire_beam_fwhm('PSW'),
            get_spire_beam_fwhm('PMW'),
            get_spire_beam_fwhm('PLW')]

    ncols = len(radave)
    increment = np.zeros(ncols)
    offsets = np.zeros(ncols)
    fit = np.zeros((2, ncols))

    for i in range(ncols):
        rc = params['rc']
        beta = params['beta']
        if not beam:
            roft = (1. + (radave[i]['midbin'] / rc)**2)**((1. - 3. * beta / 2.))
            # print(roft, beta, radave[i]['midbin'], rc)
            # plt.plot(roft)
            # plt.show()

        else:
            # create a radial profile for the beta model
            length = int(np.nanmax((radave[i]['midbin'])))
            xrad = np.arange(-1 * length, length)
            roftprime = (1. + (xrad / rc)**2)**((1. - 3. * beta) / 2.)
            pf = -1 * (2.18 * log(2) / beam[i]**2)
            psf = np.exp(pf * xrad**2)
            psf = [x / np.sum(psf) for x in psf]
            roftprimep = convolve_fft(roftprime, psf)
            roftprimep = roftprimep / np.max(np.abs(roftprimep))
            roft = np.interp(radave[i]['midbin'],xrad,roftprimep)
            plt.plot(roft)
            plt.show()
            exit()
            # roft = splrep(y=roftprimep, x=xrad, xe=radave[i]['midbin'])

        # if superplot:
            plt.plot(radave[i]['midbin'], radave[i]['fluxbin'])
            plt.show()

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
    maps,err = get_data('a0370')
    params,err = get_clus_params('a0370')
    radave = compute_rings(maps,params,30.0)
    fit = fitsz(radave, params)
