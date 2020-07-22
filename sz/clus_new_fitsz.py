################################################################################
# NAME : clus_new_fitsz
# DATE STARTED : May 25, 2020
# AUTHORS : Victoria Butler & Benjamin Vaughan
# PURPOSE :
#
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#           nsims (number of simulations/files to read in)
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################# from astropy.io import fits
import numpy as np
from math import *
import sys
sys.path.append('utilities')
sys.path.append('source_handling')
import config
from clus_get_data import clus_get_data
from FITS_tools.hcongrid import hcongrid
import matplotlib.pyplot as plt
from clus_convert_bolocam import *
from clus_get_relsz import *
from clus_get_lambdas import *
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve_fft as convolve
from scipy.optimize import curve_fit
import math

# np.random.seed(102793)

def chi_square_test(data,model,sigma):
    total_chi = []
    for i in range(len(data)):
        total_chi.append((data[i] - model[i])**2 / (sigma[i]**2))
    return np.sum(total_chi)

def fitting_func(a, x, b):
    return a*x + b

def clus_new_fitsz(maps, saveplot=0, nsim=None):

    band = ['PSW','PMW','PLW']
    fit = [0]*3

    data_file = config.CLUSDATA + 'bolocam/new/' + maps[0]['name'] + '.fits'
    base_model = fits.getdata(data_file)

    data = fits.open(data_file)
    header = data[0].header

    ylims = [[-0.01, 0.01], [0, 0.1], [0.1, 0.4]]
    for i in range(3):
        #Error Analysis -------
        error_map =  maps[i]['error']
        error_test = np.multiply(error_map, maps[i]['calfac'])
        error = error_map.flatten()

        noise_map = maps[i]['noise']

        # plt.imshow(noise_map, origin='lower')
        # plt.title('Bolocam Template_%s' % maps[i]['band'])
        # plt.colorbar()
        # plt.savefig(config.OUTPUT + 'sz_fit/fitsz_in_temp%s%s.png' % (maps[i]['band'], nsim))
        # plt.clf()
        # noise_map = np.random.normal(loc=0, scale=1.0,size=(maps[i]['signal'].shape[0],maps[i]['signal'].shape[1]))
        noise = noise_map.flatten()

        for j in range(len(noise)):
            if not math.isnan(error[j]) : #Jy/beam
                noise[j] = maps[i]['calfac']*noise[j]

        noise_map = np.reshape(noise, maps[i]['signal'].shape)

        bolo_temp = hcongrid(base_model, header, maps[i]['shead'])

        fwhm = maps[i]['widtha']
        pixscale = maps[i]['pixsize']
        retext = round(fwhm * 5.0 / pixscale)
        if retext % 2 == 0:
            retext += 1
        bmsigma = fwhm / math.sqrt(8 * math.log(2))
        beam = Gaussian2DKernel(bmsigma / pixscale, x_size=retext,y_size=retext, mode='oversample',factor=1)
        beam *= 1 / np.sum(beam.array)
        model = convolve(bolo_temp, beam, boundary='wrap')

        for j in range(maps[i]['signal'].shape[0]):
            for k in range(maps[i]['signal'].shape[1]):
                if maps[i]['mask'][j,k] == 0 :
                    model[j,k] = np.nan
                    maps[i]['srcrm'][j,k] = np.nan
                    noise_map[j,k] = np.nan
                else:
                    maps[i]['srcrm'][j,k] *= maps[i]['calfac']

        image = maps[i]['srcrm']
        #
        # plt.imshow(maps[i]['srcrm'], origin='lower')
        # plt.title('SZ without noise_ %s' % maps[i]['band'])
        # plt.colorbar()
        # plt.savefig(config.OUTPUT + 'sz_fit/sz_no_noise%s%s.png' % (maps[i]['band'], nsim))
        # plt.clf()
        #
        # plt.imshow(noise_map, origin='lower')
        # plt.title('Noise Map_%s' % maps[i]['band'])K
        # plt.colorbar()
        # plt.savefig(config.OUTPUT + 'sz_fit/noise_map%s%s.png' % (maps[i]['band'], nsim))
        # plt.clf()


        plt.imshow(image, origin='lower', clim=(-0.1, 0.1))
        plt.title('Input SZ effect %s %s' % (maps[i]['band'], nsim))
        plt.colorbar()
        plt.savefig(config.OUTPUT + 'sz_fit/fitsz_in_map%s%s.png' % (maps[i]['band'], nsim))
        plt.clf()

        plt.imshow(model, origin='lower')
        plt.title('Bolocam Template %s' % maps[i]['band'])
        plt.colorbar()
        plt.savefig(config.OUTPUT + 'sz_fit/fitsz_in_temp%s%s.png' % (maps[i]['band'], nsim))
        plt.clf()

        flat_im = image.flatten()
        flat_mod = model.flatten()

        x_data = [flat_mod[i] for i in range(len(flat_mod)) if not np.isnan(noise[i])]
        y_data = [flat_im[i] for i in range(len(flat_im)) if not np.isnan(noise[i])]
        # noise_data = [1 / n for n in noise_map.flatten() if not np.isnan(n)] #polyfit wants weights as 1 / sigma
        sigma_data = [n for n in noise if not np.isnan(n)]


        print(len(x_data), len(y_data), len(sigma_data))
        # z, cov = np.polyfit(x_data,y_data,1, cov=True),# w = noise_data)

        params, cov = curve_fit(fitting_func, x_data, y_data, sigma=sigma_data)
        e = np.sqrt(np.diag(cov))
        # p = np.poly1d(z)
        intercept = params[1]
        intercept_error = e[1]
        slope = params[0]
        slope_error = e[0]
        y_fit = [slope * x + intercept for x in x_data]
        fit[i] = params
        chi_square = chi_square_test(y_data ,y_fit, sigma_data)
        red_chi_square = chi_square / len(y_data)


        if saveplot :
            input_DI = np.load(config.OUTPUT + 'add_sz/input_dI_%s.npy' %(band[i]))
            # plt.ylim(ylims[i])
            plt.plot(x_data, y_fit, c='red', label='dI = %.4E +/- %.3E \n offset = %.1E +/- %.1E \n chi square: %.3E \n red chi square %.3E' % (slope, slope_error, intercept, intercept_error, chi_square, red_chi_square))
            plt.scatter(x_data,y_data, alpha=0.5, label='Image Data: Input DI : %.4E' % input_DI)
            in_y = [input_DI * x  for x in x_data]
            # plt.plot(x_data, in_y, label='input dI = %.4E' % input_DI, c='green')

            sigma = (slope - input_DI) / slope_error
            print('uncertainty in fits', slope_error, intercept_error)
            print( 'Slope is %.2E sigma from the expected value' % sigma )

            #these next few lines is for testing only

            plt.legend()
            plt.xlabel('Model Flux [Unitless]')
            plt.ylabel('Image Flux [MJy/Sr]')
            plt.title('SZ Fit for %s %s' % (maps[i]['band'], nsim))
            # plt.xlim(-0.005,0.005)
            plt.savefig(config.OUTPUT + 'sz_fit/rxj1347_%s_%s.png' % (maps[i]['band'], nsim)) #rxj1347 placeholder for clusname.
            plt.clf()

    return fit

if __name__ == '__main__' :
    maps,err = clus_get_data('rxj1347','0',sgen=3)
    sim = 0

    lense_file = config.HOME + 'Lensing/lense_template_' + band[i] + '.fits'
    lense_model = fits.getdata(lense_file)

    # file = config.OUTPUT + 'pcat_residuals/rxj1347_mask_resid_%s_%s.fits'%(band[i],sim)
    file = config.OUTPUT + 'pcat_residuals/rxj1347_resid_%s_%s.fits'%(band[i],sim)
    image = fits.getdata(file) - lense_model

    clus_new_fitsz(image)
