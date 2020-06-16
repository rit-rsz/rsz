##########################################################################
# NAME : clus_popmap.py
# DATE STARTED : October 3, 2019
# AUTHORS : Victoria Butler
# PURPOSE : This program populates simulated maps with sources lensed by lenstool.
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :

# OUTPUTS :
# REVISION HISTORY :

##########################################################################
from astropy.wcs.utils import skycoord_to_pixel, pixel_to_skycoord
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import sys, os
sys.path.append('../utilities')
import config
from gaussian import makeGaussian
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import math
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve_fft as convolve
sys.path.append('../multiband_pcat')
from image_eval import psf_poly_fit, image_model_eval
import pdb

def clus_popmap(ltfile,maps,map_size,band,name,pixsize,fwhm,loz=None,superplot=0,savemaps=0,genbethermin=None):

    # read in the image.dat file
    ra = []
    dec = []
    mag = []
    with open(ltfile,'r') as f :
        data = f.readlines()
        for i in range(len(data)) :
            val = data[i].strip().split(' ')
            if i == 0 :
                racent = float(val[-2])
                deccent = float(val[-1])
            else :
                ra.append(float(val[1]))
                dec.append(float(val[2]))
                mag.append(float(val[-1]))
    f.close()

    # add low redshift sources back into map
    print(len(loz['x']))
    ra.extend(loz['x'])
    dec.extend(loz['y'])
    flux = [10.0**(-k/2.5) for k in mag]
    flux.extend(loz['f'])

    refx = float(racent)
    refy = float(deccent)

    # #convert ra/dec to degrees
    '''
    Do a sign check on central RA/DEC coordinates
    '''
    ra = [((-x / 3600.0) + refx) for x in ra]
    dec = [((y / 3600.0) + refy) for y in dec]

    if superplot:
        plt.scatter(ra,dec,s=2,c=flux)
        plt.colorbar()
        plt.title('Clus Popmap: Pixels')
        plt.show()
        plt.clf()

    header = maps['shead']
    wcs = WCS(header)
    coords = SkyCoord(ra=ra*u.deg,dec=dec*u.deg)
    x,y = skycoord_to_pixel(coords, wcs)

    truthtable = {'x': x, 'y': y,'flux': flux}

    if superplot:
        plt.scatter(x,y,s=2,c=flux)
        plt.colorbar()
        plt.title('Clus Popmap: Skycoord')
        plt.show()
        plt.clf()

    ''' PCAT LION METHOD'''
    orig_length = len(x)
    mapy = maps['signal'].shape[0]
    mapx = maps['signal'].shape[1]
    x = np.array(x,dtype=np.float32)
    y = np.array(y,dtype=np.float32)
    flux = np.array(flux,dtype=np.float32)
    pixel_per_beam = 2*np.pi*((3)/2.355)**2
    psf, cf, nc, nbin = get_gaussian_psf_template(fwhm,pixel_fwhm=3.) # assumes pixel fwhm is 3 pixels in each band
    sim_map = image_model_eval(y, x, pixel_per_beam*nc*flux, 0.0, (mapx, mapy), int(nc), cf)
    '''#############################################################################################################'''

    # if superplot:
    # plt.imshow(sim_map,extent=(0,300,0,300),clim=[0.0,0.15],origin=0)
    plt.imshow(sim_map,clim=[0,0.01],origin=0)
    plt.colorbar()
    plt.title('Clus Popmap: Lensed map')
    plt.savefig('lensedmap_' + name + '_' + band + '.png')
    plt.clf()

    if savemaps:
        hdx = fits.PrimaryHDU(maps['signal'],maps['shead'])
        sz = fits.PrimaryHDU(sim_map,hdx.header)
        # sz.writeto(config.SIMBOX + 'lensedmap_' + name + '_' + band + '.fits',overwrite=True)
        sz.writeto(config.SIM + 'lensedmap_' + name + '_' + band + '.fits',overwrite=True)

    return sim_map, truthtable

def get_gaussian_psf_template(fwhm,pixel_fwhm=3., nbin=5):
    nc = nbin**2
    psfnew = Gaussian2DKernel(pixel_fwhm/2.355*nbin, x_size=125, y_size=125).array.astype(np.float32)
    # psfnew2 = psfnew / np.max(psfnew)  * nc
    cf = psf_poly_fit(psfnew, nbin=nbin)
    return psfnew, cf, nc, nbin
