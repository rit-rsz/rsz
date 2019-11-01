################################################################################
# NAME : clus_fitsz.py
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
from clus_get_clusparams import *
from clus_compute_rings import *
from clus_get_data import *
from scipy.stats import linregress
from astropy.convolution import convolve_fft
from scipy.interpolate import splrep
from scipy.optimize import curve_fit
from gaussian import *
from get_spire_beam_fwhm import *

def fitting_func(a, x, b):
    return a*x + b

def clus_fitsz(radave, params, beam=None, maxlim=3600, minlim=0, noweight=1, superplot=1, verbose=1, nsim=0, saveplot=1):

    # init params
    ncols = len(radave)
    increment = np.zeros(ncols)
    offsets = np.zeros(ncols)
    fit = [0]*3
    rc = params['rc']
    beta = params['beta']
    for i in range(ncols):

        if not beam:
            roft = (1. + (radave[i]['midbin'] / rc)**2)**((1. - 3. * beta / 2.))

        else:
            # create a radial profile for the beta model
            length = int(np.nanmax((radave[i]['midbin'])))
            xrad = np.arange(-1 * length, length)
            roftprime = (1. + (xrad / rc)**2)**((1. - 3. * beta) / 2.)
            #create psf for beam
            sigma = beam[i] / sqrt(8 * log(2))
            pf = -1 * (1 / (2 * sigma**2))
            psf = np.exp(pf * xrad**2)
            psf = [x / np.sum(psf) for x in psf]
            #convolve beam's psf with r(t)
            roftprimep = convolve_fft(roftprime, psf)
            roftprimep = roftprimep / np.max(np.abs(roftprimep))
            roft = np.interp(radave[i]['midbin'],xrad,roftprimep)

        if superplot or saveplot:
            plt.plot(radave[i]['midbin'], radave[i]['fluxbin'])
            plt.ylabel('Radial Average (MJy/sr)')
            plt.xlabel('Radius (arcsec)')
            plt.title(params['clusname'] + '  ' + radave[i]['band'] + '  RC:' + str(rc) + '  Beta:' + str(beta))
            plt.xlim((0,600))
            # plt.ylim((0.0,0.3))
            plt.errorbar(radave[i]['midbin'],radave[i]['fluxbin'],yerr=radave[i]['errbin'])
            plt.scatter(radave[i]['midbin'], radave[i]['fluxbin'],marker='o')
            if superplot:
                plt.show()
            elif saveplot:
                if nsim != 0:
                    filename = config.HOME + 'outputs/IB_fits/ra_v_r_' + params['clusname'] + '_' + radave[i]['band'] + '_' + str(nsim) + '.pdf'
                else:
                    filename = config.HOME + 'outputs/IB_fits/ra_v_r_' + params['clusname'] + '_' + radave[i]['band'] + '_' + '.pdf'
                plt.savefig(filename, format='pdf')

        # ignore values above radius of 600
        index = []
        for k in range(len(radave[i]['midbin'])):
            if radave[i]['midbin'][k] >= 600:
                index.append(k)
        radave[i]['midbin'] = np.delete(radave[i]['midbin'], k)
        radave[i]['fluxbin'] = np.delete(radave[i]['fluxbin'], k)
        roft = np.delete(roft, k)

        # remove nans from the dataset
        for k in range(len(radave[i]['midbin'])):
            if np.isnan(radave[i]['midbin'][k]) or np.isnan(radave[i]['fluxbin'][k]):
                radave[i]['midbin'][k] = 0
                radave[i]['fluxbin'][k] = 0
                roft[k] = 0

        #create a linear fit for r(t) vs fluxbin
        z, cov = curve_fit(fitting_func, roft, radave[i]['fluxbin'])
        intercept = z[1]
        slope = z[0]
        fit[i] = z
        line = slope * roft + intercept

        if superplot or saveplot:
            plt.plot(roft, radave[i]['fluxbin'],label='Data')
            plt.plot(roft, line,label='Best Linear Fit')
            plt.xlabel(r'R($\theta$)')
            plt.ylabel('Radial Average (MJy/sr)')
            plt.legend()

            # ax[1].set_xlabel(r"R($\theta^\prime$)")
            # ax[1].set_ylabel('Radial Average (MJy/sr)')
            # ax[1].plot(xrad, roftprimep,label='Data')
            # ax[1].plot(xrad, p(roftprimep),label='Best Linear Fit')
            # ax[1].set_xlabel(r"R($\theta^\prime$)")
            # ax[1].set_ylabel('Radial Average (MJy/sr)')
            # ax[1].legend()

            plt.title(params['clusname'] + '  ' + radave[i]['band'] + '  Slope: %.4f  Intercept: %.4f' %(slope,intercept))
            # fig.tight_layout()
            if superplot:
                plt.show()
            if saveplot:
                if nsim != 0:
                    filename = config.HOME + 'outputs/IB_fits/dI_fit_' + params['clusname'] + '_' + radave[i]['band'] + '_' + str(nsim) + '.pdf'
                else:
                    filename = config.HOME + 'outputs/IB_fits/dI_fit_' + params['clusname'] + '_' + radave[i]['band'] + '_' + '.pdf'
                plt.savefig(filename, format='pdf')

    return fit

if __name__ == '__main__':
    maps,err = clus_get_data('a0370')
    params,err = clus_get_clusparams('a0370')
    radave = clus_compute_rings(maps,params,30.0)
    fit = clus_fitsz(maps,radave, params)
