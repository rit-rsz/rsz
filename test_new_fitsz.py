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

band = ['PSW','PMW','PLW']
maps,err = clus_get_data('rxj1347','0',sgen=3)
sim = 0

for i in range(3):
    lense_file = config.HOME + 'Lensing/lense_template_' + band[i] + '.fits'
    lense_model = fits.getdata(lense_file)

    # file = config.OUTPUT + 'pcat_residuals/rxj1347_mask_resid_%s_%s.fits'%(band[i],sim)
    file = config.OUTPUT + 'pcat_residuals/rxj1347_resid_%s_%s.fits'%(band[i],sim)
    image = fits.getdata(file) - lense_model
    # image = fits.getdata(file)

    file2 = config.HOME + 'sz/bolocam_spire_%s.fits'%(band[i])
    model = fits.getdata(file2)


    for j in range(maps[i]['signal'].shape[0]):
        for k in range(maps[i]['signal'].shape[1]):
            if maps[i]['mask'][j,k] == 0 :
                model[j,k] = np.nan
                image[j,k] = np.nan

    if i == 0 :
        plw_image = hcongrid(image,maps[0]['shead'],maps[2]['shead'])
        plw_model = hcongrid(model,maps[0]['shead'],maps[2]['shead'])
    elif i == 1 :
        plw_image = hcongrid(image,maps[1]['shead'],maps[2]['shead'])
        plw_model = hcongrid(model,maps[1]['shead'],maps[2]['shead'])

    else :
        plw_image = image
        plw_model = model


    x_data = [x*1e3 for x in plw_model.flatten() if not np.isnan(x)]
    y_data = [y for y in plw_image.flatten() if not np.isnan(y)]
    zero_mask = x_data != 0
    x_data[zero_mask]
    y_data[zero_mask]


    z = np.polyfit(x_data,y_data,1)
    print(z)
    p = np.poly1d(z)
    intercept = z[1]
    slope = z[0]
    plt.plot(x_data,p(x_data))
    plt.scatter(x_data,y_data)
    # plt.xlim(-0.005,0.005)
    plt.ylim(-0.05,0.05)
    plt.savefig('fitsz_test' + band[i] + '.png')
    plt.clf()
