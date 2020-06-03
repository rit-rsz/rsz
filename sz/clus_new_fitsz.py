from astropy.io import fits
import numpy as np
from math import *
import sys
sys.path.append('utilities')
sys.path.append('source_handling')
import config
from clus_get_data import clus_get_data
from FITS_tools.hcongrid import hcongrid
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def fitting_func(a, x, b):
    return a*x + b

def clus_new_fitsz(maps, saveplot=0):

    band = ['PSW','PMW','PLW']
    fit = [0]*3

    for i in range(3):
        bol_temp = config.OUTPUT + 'bolocam_spire_%s.fits'%(band[i])
        model = fits.getdata(bol_temp)

        # for j in range(maps[i]['signal'].shape[0]):
        #     for k in range(maps[i]['signal'].shape[1]):
        #         if maps[i]['mask'][j,k] == 0 :
        #             model[j,k] = np.nan
        #             maps[i]['srcrm'][j,k] = np.nan

        mask = maps[i]['mask'] == 0
        maps[i]['srcrm'][mask] = np.nan
        model[mask] = np.nan

        if i == 0 :
            image = hcongrid(maps[i]['srcrm'],maps[0]['shead'],maps[2]['shead'])
            model = hcongrid(model,maps[0]['shead'],maps[2]['shead'])
        elif i == 1 :
            image = hcongrid(maps[i]['srcrm'],maps[1]['shead'],maps[2]['shead'])
            model = hcongrid(model,maps[1]['shead'],maps[2]['shead'])

        else :
            image = maps[i]['srcrm']
            model = model


        x_data = [x*1e3 for x in model.flatten() if not np.isnan(x)]
        y_data = [y for y in image.flatten() if not np.isnan(y)]
        zero_mask = x_data != 0
        x_data[zero_mask]
        y_data[zero_mask]


        z = np.polyfit(x_data,y_data,1)
        print(z)
        p = np.poly1d(z)
        intercept = z[1]
        slope = z[0]
        fit[i] = z

        if saveplot :
            plt.plot(x_data,p(x_data))
            plt.scatter(x_data,y_data)
            # plt.xlim(-0.005,0.005)
            plt.ylim(-0.05,0.05)
            plt.savefig('fitsz_test' + band[i] + '.png')
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
