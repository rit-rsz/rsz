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
from gaussian import makeGaussian
import matplotlib.pyplot as plt
from astropy.io import fits
from FITS_tools.hcongrid import hcongrid
import numpy as np
import math
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve_fft as convolve
sys.path.append('../multiband_pcat/multiband_pcat')
from image_eval import psf_poly_fit, image_model_eval

def clus_popmap(ltfile,maps,map_size,band,name,pixsize,fwhm,loz=None,superplot=0,savemaps=0):

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

    # if len(loz) == 8 :
    #     ra = [ra,loz['x']]
    #     dec = [dec,loz['y']]

    refx = float(racent)
    refy = float(deccent)

    # convert ra/dec to degrees
    ra = [((-x / 3600.0) + refx) for x in ra]
    dec = [((y / 3600.0) + refy) for y in dec]
    flux = [10.0**(-k/2.5) for k in mag]

    # if len(loz) == 8 :
    #     flux = [flux,loz['f']]

    if superplot:
        plt.scatter(ra,dec,s=2,c=flux)
        plt.colorbar()
        plt.title('Clus Popmap: Pixels')
        plt.show()

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

    # make the output map to have noise added later
    # xf = np.floor(x)
    # yf = np.floor(y)
    # cmap = np.zeros((map_size,map_size),dtype=np.float32)
    # nx, ny = cmap.shape
    # np.place(xf, xf > nx-1, nx-1)
    # np.place(yf, yf > ny-1, ny-1)
    # for cx, cy, cf in zip(xf, yf, flux):
    #     cmap[int(cy), int(cx)] += cf  # Note transpose
    #
    # beam = get_gauss_beam(fwhm,pixsize,band,oversamp=5)
    # sim_map = convolve(cmap, beam, boundary='wrap')

    ''' PCAT LION METHOD'''
    orig_length = len(x)
    position_mask = np.logical_and(np.logical_and(x > 0, x < map_size), np.logical_and(y > 0, y < map_size))
    x = np.array(x,dtype=np.float32)[position_mask]
    y = np.array(y,dtype=np.float32)[position_mask]
    flux = np.array(flux,dtype=np.float32)[position_mask]

    psf, cf, nc, nbin = get_gaussian_psf_template(fwhm,pixel_fwhm=3.) # assumes pixel fwhm is 3 pixels in each band
    sim_map = image_model_eval(x, y, nc*flux, 0.0, (int(map_size), int(map_size)), int(nc), cf)
    '''#############################################################################################################'''

    if superplot:
        plt.imshow(sim_map,extent=(0,300,0,300),clim=[0.0,0.15],origin=0)
        plt.colorbar()
        plt.title('Clus Popmap: Lensed map')
        plt.show()

    if savemaps:
        hdx = fits.PrimaryHDU(maps['signal'],maps['shead'])
        sz = fits.PrimaryHDU(sim_map,hdx.header)
        if os.path.isfile(config.SIMBOX + 'lensedmap_' + name + '_' + band + '.fits'):
            os.remove(config.SIMBOX + 'lensedmap_' + name + '_' + band + '.fits')
        sz.writeto(config.SIMBOX + 'lensedmap_' + name + '_' + band + '.fits')

    return sim_map, truthtable

def get_gauss_beam(fwhm, pixscale, band, nfwhm=5.0, oversamp=1):
    retext = round(fwhm * nfwhm / pixscale)
    if retext % 2 == 0:
        retext += 1

    bmsigma = fwhm / math.sqrt(8 * math.log(2))

    beam = Gaussian2DKernel(bmsigma / pixscale, x_size=retext,
                            y_size=retext, mode='oversample',
                            factor=oversamp)
    beam *= 1.0 / beam.array.max()
    return beam

def get_gaussian_psf_template(fwhm,pixel_fwhm=3., nbin=5):
    nc = nbin**2
    psfnew = Gaussian2DKernel(pixel_fwhm/2.355*nbin, x_size=125, y_size=125).array.astype(np.float32)
    psfnew2 = psfnew / np.max(psfnew) * nc
    cf = psf_poly_fit(psfnew, nbin=nbin)
    return psfnew2, cf, nc, nbin

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
    pixsize = [6.0, 8.33333, 12.0]
    maps, err = clus_get_data(clusname=clusname, resolution=resolution, bolocam=bolocam, verbose=verbose)
    band = maps[0]['band']
    ltfile = config.SIMBOX + 'a0370' + '_image_' + band + '.dat'
    # ltfile = '/data/zemcov/clusters/sandbox/a0370_image_PSW.dat'
    clus_popmap(ltfile,maps[0],band,'a0370',pixsize[0],superplot=1,savemaps=1)
