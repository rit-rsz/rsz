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
from scipy.optimize import curve_fit
from gaussian import *
from get_spire_beam_fwhm import *

def fitting_func(a, x, b):
    return a*x + b

def clus_fitsz(radave, params, beam, sgen=None, maxlim=3600, minlim=0, noweight=1, superplot=1, verbose=1, nsim=0, saveplot=1):

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
            plt.title('Clus Fitsz : ' + params['clusname'] + '  ' + radave[i]['band'] + '  RC:' + str(rc) + '  Beta:' + str(beta))
            plt.xlim((0,600))
            plt.errorbar(radave[i]['midbin'],radave[i]['fluxbin'],yerr=radave[i]['errbin'])
            plt.scatter(radave[i]['midbin'], radave[i]['fluxbin'],marker='o')
            if superplot:
                plt.show()
            elif saveplot:
                if sgen != None:
                    filename = config.OUTPUT + 'radial_averages/' + params['clusname'] + '_ra_v_r_' + radave[i]['band'] + '_' + str(nsim) + '.png'
                else:
                    filename = config.OUTPUT + 'radial_averages/' + params['clusname'] + '_ra_v_r_' + radave[i]['band'] + '_real.png'

                plt.savefig(filename)
                plt.clf()

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
        # np.place(roft, roft<0.1, 0.0)
        # roft = [x for x in roft if x != 0]
        # radave[i]['fluxbin'] = [y for y in radave[i]['fluxbin'] if y != 0][:len(roft)]
        radave[i]['fluxbin'] = radave[i]['fluxbin'][:len(roft)]
        print(radave[i]['fluxbin'])
        # z, cov = curve_fit(fitting_func, roft, radave[i]['fluxbin'],sigma=radave[i]['errbin'])
        z, cov = curve_fit(fitting_func, roft, radave[i]['fluxbin'])
        intercept = z[1]
        slope = z[0]
        fit[i] = z
        line = [(slope*x) + intercept for x in roft]

        if superplot or saveplot:
            plt.scatter(roft, radave[i]['fluxbin'],label='Data')
            plt.plot(roft, line,label='Best Linear Fit')
            plt.xlabel(r'R($\theta$)')
            plt.ylabel('Radial Average (MJy/sr)')
            # plt.ylim((-0.1,0.1))
            plt.legend()
            plt.title('Clus Fitsz : ' + params['clusname'] + '  ' + radave[i]['band'] + '  Slope: %.4f  Intercept: %.4f' %(slope,intercept))
            if superplot:
                plt.show()
            if saveplot:
                if sgen != None:
                    filename = config.OUTPUT + 'IB_fits/' + params['clusname'] + '_dI_fit_' + radave[i]['band'] + '_' + str(nsim) + '.png'
                else:
                    filename = config.OUTPUT + 'IB_fits/' + params['clusname'] + '_dI_fit_' + radave[i]['band'] + '_' + 'real.png'

                plt.savefig(filename)
                plt.clf()

    return fit
