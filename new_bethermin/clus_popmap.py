
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
from astropy.wcs.utils import skycoord_to_pixel
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import sys, os
sys.path.append('../utilities')
import config
from gaussian import makeGaussian, padGaussian
from astropy.convolution import Gaussian2DKernel
import matplotlib.pyplot as plt
from astropy.io import fits
from FITS_tools.hcongrid import hcongrid
import numpy as np
sys.path.append('../multiband_pcat')
from image_eval import psf_poly_fit, image_model_eval
from astropy.convolution import convolve_fft as convolve
from scipy.ndimage import shift
import time
import math

def clus_popmap(ltfile,maps,band,name,pixsize,fwhm,loz=None,superplot=0,savemaps=0):

    # read in the image.dat file
    ra = []
    dec = []
    mag = []
    with open(ltfile,'r') as f :
        data = f.readlines()
        print('length of popmap data:',len(data))
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

    refx = float(racent)
    refy = float(deccent)
    print('popmap center vals:',refx,refy)
    # convert ra/dec to degrees
    ra = [((-x / 3600.0) + refx) for x in ra]
    dec = [((y / 3600.0) + refy) for y in dec]
    flux = [10.0**(-z/2.5) for z in mag]

    #stick our cropped sources back into lensed array
    # os.system('say "program finished"')
    # print('\a')
    #
    # ra.extend([((-1.0*j / 3600.0) + refx) for j in loz[0]])
    # dec.extend([((k / 3600.0) + refy) for k in loz[1]])
    # flux.extend(loz[3])
    # flux.extend([10.0**(-l/2.5) for l in loz[3]])

    # if superplot:
    plt.scatter(ra,dec,s=2,c=flux)
    plt.colorbar()
    plt.title('Start of Clus Popmap')
    plt.savefig('popmap_%s' %(band))
    plt.clf()

    header = maps['shead']
    wcs = WCS(header)
    coords = SkyCoord(ra=ra*u.deg,dec=dec*u.deg)
    x,y = skycoord_to_pixel(coords, wcs)

    # if superplot:
    plt.scatter(x,y,s=2,c=flux)
    plt.colorbar()
    plt.title('clus_popmap: after skycoord_to_pixel')
    plt.savefig('popmap_end_%s' %(band))
    plt.clf()

    x_size = maps['signal'].shape[0]
    y_size = maps['signal'].shape[1]
    pixscale = [6.0, 8.33333, 12.0]
    outmap = np.zeros((x_size,y_size),dtype=np.float32)
    # cmap = np.zeros((300,300),dtype=np.float32)

    # if band == 'PSW':
    #     # xf = np.floor(xpos)
    #     # yf = np.floor(ypos)
    #     xf = np.floor(x)
    #     yf = np.floor(y)
    #     nx, ny = cmap.shape
    #     np.place(xf, xf > nx-1, nx-1)
    #     np.place(yf, yf > ny-1, ny-1)
    #     for cx, cy, cf in zip(xf, yf, flux):
    #         cmap[int(cy), int(cx)] += cf  # Note transpose
    #
    # # Other bands, with pixel scale adjustment
    # else :
    #     if band == 'PMW':
    #         posrescale = pixscale[0] / pixscale[1]
    #     if band == 'PLW':
    #         posrescale = pixscale[0] / pixscale[2]
    #
    #     xf = np.floor(posrescale * x)
    #     yf = np.floor(posrescale * y)
    #     nx, ny = cmap.shape
    #     np.place(xf, xf > nx-1, nx-1)
    #     np.place(yf, yf > ny-1, ny-1)
    #     for cx, cy, cf in zip(xf, yf, flux):
    #         cmap[int(cy), int(cx)] += cf  # Note transpose


    num = 0
    fig, axs = plt.subplots(1, 2)

    # print(datetime.now(), 'Time before manually multiplying through gaussians')
    for i in range(len(flux)):
        if x[i] > 0 and x[i] < y_size and y[i] > 0 and y[i] < x_size:
            sigma,kern = makeGaussian(y_size,x_size, fwhm = fwhm/pixsize, center=(x[i],y[i]))
            norm = np.sum(kern)**2
            psf = kern * flux[i] * norm
            outmap = outmap + psf
        else :
            num += 1

    print('number of sources outside map: ',num)

    # beam = get_gauss_beam(fwhm,pixsize,band)
    # plt.imshow(cmap)
    # plt.savefig('cmap_%s.png' %(band))
    # sim_map = convolve(cmap, beam, boundary='wrap')

    # orig_length = len(x)
    # print(len(x),len(y),len(flux))
    # print(len(loz[0]),len(loz[1]),len(loz[3]))
    # position_mask = np.logical_and(np.logical_and(x > 0, x < y_size), np.logical_and(y > 0, y < x_size))
    # print('position mask: ',len(position_mask))
    # x = np.array(x,dtype=np.float32)[position_mask]
    # y = np.array(y,dtype=np.float32)[position_mask]
    # flux = np.array(flux,dtype=np.float32)[position_mask]
    #
    # psf, cf, nc, nbin = get_gaussian_psf_template(fwhm,pixel_fwhm=3.) # assumes pixel fwhm is 3 pixels in each band
    # outmap = image_model_eval(x, y, nc*flux, 0.0, (int(y_size), int(x_size)), int(nc), cf)

    # outmap = clus_image_model(fwhm/pixsize,input_outmap,x,y,flux)

    if superplot:
        plt.imshow(outmap,extent=(0,300,0,300),clim=[0.0,0.15],origin=0)
        plt.colorbar()
        plt.title('clus_popmap: Lensed map')
        plt.show()

    if savemaps == 2:
        hdx = fits.PrimaryHDU(maps['signal'],maps['shead'])
        sz = fits.PrimaryHDU(outmap,hdx.header)
        sz.writeto(config.SIMBOX + 'lensedmap_' + name + '_' + band + '.fits',overwrite=True)

    # return sim_map,beam
    return outmap

def get_gaussian_psf_template(fwhm,pixel_fwhm=3., nbin=5):
    nc = nbin**2
    psfnew = Gaussian2DKernel(pixel_fwhm/2.355*nbin, x_size=125, y_size=125).array.astype(np.float32)
    psfnew2 = psfnew / np.max(psfnew) * nc
    cf = psf_poly_fit(psfnew, nbin=nbin)
    return psfnew2, cf, nc, nbin

def get_gauss_beam(fwhm, pixscale, band, nfwhm=5.0, oversamp=1):
    print('get_gauss_beam: fwhm',fwhm)
    print('get_gauss_beam: pixscale',pixscale)
    retext = round(fwhm * nfwhm / pixscale)
    if retext % 2 == 0:
        retext += 1

    bmsigma = fwhm / math.sqrt(8 * math.log(2))

    beam = Gaussian2DKernel(bmsigma / pixscale, x_size=retext,
                            y_size=retext, mode='oversample',
                            factor=oversamp)
    beam *= 1.0 / beam.array.max()
    plt.imshow(beam)
    plt.savefig('beam_%s.png' %(band))
    return beam

if __name__ == '__main__':
    import sys
    sys.path.append('../utilities')
    sys.path.append('../source_handling')
    from clus_get_data import clus_get_data
    import config
    clusname = 'a0370'
    resolution = 'nr'
    bolocam = None
    verbose = 1
    fwhm = [17.6, 23.9, 35.2]
    pixsize = [6.0, 8.33333, 12.0]
    maps, err = clus_get_data(clusname=clusname, resolution=resolution, bolocam=bolocam, verbose=verbose)
    band = maps[0]['band']
    # ltfile = config.SIMBOX + 'a0370' + '_image_' + band + '.dat'
    ltfile = config.SIMBOX + 'a0370_image_PSW_test.dat'
    clus_popmap(ltfile,maps[0],band,'a0370',pixsize[0],fwhm[0],superplot=0,savemaps=1)
