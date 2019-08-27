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
from scipy.optimize import curve_fit
from gaussian import *
from get_spire_beam_fwhm import *

def fitsz(radave, params, beam=None, maxlim=3600, minlim=0, noweight=1, superplot=1, verbose=1):

    # init params
    # beam = [get_spire_beam_fwhm('PSW'),
    #         get_spire_beam_fwhm('PMW'),
    #         get_spire_beam_fwhm('PLW')]

    ncols = len(radave)
    increment = np.zeros(ncols)
    offsets = np.zeros(ncols)
    fit = [0]*3

    for i in range(ncols):
        rc = params['rc']
        beta = params['beta']
        if not beam:
            roft = (1. + (radave[i]['midbin'] / rc)**2)**((1. - 3. * beta / 2.))

        else:
            # create a radial profile for the beta model
            length = int(np.nanmax((radave[i]['midbin'])))
            xrad = np.arange(-1 * length, length)
            roftprime = (1. + (xrad / rc)**2)**((1. - 3. * beta) / 2.)
            #create psf for beam
            pf = -1 * (2.18 * log(2) / beam[i]**2)
            psf = np.exp(pf * xrad**2)
            psf = [x / np.sum(psf) for x in psf]
            #convolve beam's psf with r(t)
            roftprimep = convolve_fft(roftprime, psf)
            roftprimep = roftprimep / np.max(np.abs(roftprimep))
            roft = np.interp(radave[i]['midbin'],xrad,roftprimep)

        if superplot:
            plt.plot(radave[i]['midbin'], radave[i]['fluxbin'])
            plt.ylabel('Radial Average (MJy/sr)')
            plt.title(params['clusname'] + '  ' + radave[i]['band'] + '  RC:' + str(rc) + '  Beta:' + str(beta))
            plt.xlim((0,600))
            plt.ylim((-2,2))
            plt.errorbar(radave[i]['midbin'],radave[i]['fluxbin'],yerr=radave[i]['errbin'])
            plt.scatter(radave[i]['midbin'], radave[i]['fluxbin'],marker='o')
            plt.show()

        # ignore values above radius of 600
        # radave[i]['midbin'] = radave[i]['midbin'][15:]
        # radave[i]['fluxbin'] = radave[i]['fluxbin'][15:]



        # remove nans from the dataset
        for k in range(len(radave[i]['midbin'])):
            if np.isnan(radave[i]['midbin'][k]) :
                radave[i]['midbin'][k] = 0
                radave[i]['fluxbin'][k] = 0
                roft[k] = 0
                roftprimep[k] = 0

        # slope, intercept, r_value, p_value, std_err = linregress(roft, radave[i]['fluxbin'])
        #create a linear fit for r(t) vs fluxbin
        z = np.polyfit(roft,radave[i]['fluxbin'],1)
        p = np.poly1d(z)
        intercept = z[1]
        slope = z[0]
        fit[i] = z
        print(z)

        if superplot:
            fig, ax = plt.subplots(nrows=1 , ncols=2)
            ax[0].plot(roft, radave[i]['fluxbin'],label='Data')
            ax[0].plot(roft, p(roft),label='Best Linear Fit')
            ax[0].set_xlabel(r'R($\theta$)')
            ax[0].set_ylabel('Radial Average (MJy/sr)')
            ax[0].legend()

            ax[1].set_xlabel(r"R($\theta^\prime$)")
            ax[1].set_ylabel('Radial Average (MJy/sr)')
            ax[1].plot(xrad, roftprimep,label='Data')
            ax[1].plot(xrad, p(roftprimep),label='Best Linear Fit')
            ax[1].set_xlabel(r"R($\theta^\prime$)")
            ax[1].set_ylabel('Radial Average (MJy/sr)')
            ax[1].legend()

            fig.suptitle(params['clusname'] + '  ' + radave[i]['band'] + '  Slope: %.4f  Intercept: %.4f' %(slope,intercept))
            # fig.tight_layout()
            plt.show()

    return fit

if __name__ == '__main__':
    maps,err = get_data('a0370')
    params,err = get_clus_params('a0370')
    radave = compute_rings(maps,params,30.0)
    fit = fitsz(maps,radave, params)
