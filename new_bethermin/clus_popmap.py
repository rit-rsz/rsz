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

    print('mags:',mag[0:10])
    # if len(loz) == 8 :
    #     ra = [ra,loz['x']]
    #     dec = [dec,loz['y']]

    # refx = maps['shead']['CRVAL1']
    # refy = maps['shead']['CRVAL2']
    refx = float(racent)
    refy = float(deccent)

    # convert ra/dec to degrees
    ra = [((-x / 3600.0) + refx) for x in ra]
    dec = [((y / 3600.0) + refy) for y in dec]
    flux = [10.0**(-k/2.5) for k in mag]

    # if len(loz) == 8 :
    #     flux = [flux,loz['f']]

    # if superplot:
    plt.scatter(ra,dec,s=2,c=flux)
    plt.colorbar()
    plt.xlim((39.25,40.25))
    plt.ylim((-2.0,-1.1))
    plt.title('Start of Clus Popmap')
    plt.savefig('popmap_beg_%s' %(band))
    plt.clf()

    header = maps['shead']
    wcs = WCS(header)
    print('len of data array before',len(ra),len(dec))
    coords = SkyCoord(ra=ra*u.deg,dec=dec*u.deg)
    x,y = skycoord_to_pixel(coords, wcs)
    print('len of data array after',len(x),len(y))

    # if superplot:
    plt.scatter(x,y,s=2,c=flux)
    plt.colorbar()
    plt.title('clus_popmap: after skycoord_to_pixel')
    plt.savefig('popmap_end_%s' %(band))
    plt.clf()

    x_size = maps['signal'].shape[0]
    y_size = maps['signal'].shape[1]
    # x_size = 300
    # y_size = 300
    # outmap = np.zeros((x_size,y_size),dtype=np.float32)
    # num = 0
    # for i in range(len(flux)):
    #     if x[i] > 0 and x[i] < y_size and y[i] > 0 and y[i] < x_size:
    #         sigma,kern = makeGaussian(y_size,x_size, fwhm = fwhm/pixsize, center=(x[i],y[i]))
    #         norm = kern / np.max(kern)
    #         psf = norm * flux[i]
    #         outmap = outmap + psf
    #         print(i,len(flux))
    #     else :
    #         num += 1
    # print('number of sources outside map: ',num)

    pixscale = [6.0, 8.33333, 12.0]
    cmap = np.zeros((300,300),dtype=np.float32)

    if band == 'PSW':
        # xf = np.floor(xpos)
        # yf = np.floor(ypos)
        xf = np.floor(x)
        yf = np.floor(y)
        nx, ny = cmap.shape
        np.place(xf, xf > nx-1, nx-1)
        np.place(yf, yf > ny-1, ny-1)
        for cx, cy, cf in zip(xf, yf, flux):
            cmap[int(cy), int(cx)] += cf  # Note transpose

    # Other bands, with pixel scale adjustment
    else :
        if band == 'PMW':
            posrescale = pixscale[0] / pixscale[1]
        if band == 'PLW':
            posrescale = pixscale[0] / pixscale[2]

        xf = np.floor(posrescale * x)
        yf = np.floor(posrescale * y)
        nx, ny = cmap.shape
        np.place(xf, xf > nx-1, nx-1)
        np.place(yf, yf > ny-1, ny-1)
        for cx, cy, cf in zip(xf, yf, flux):
            cmap[int(cy), int(cx)] += cf  # Note transpose
        del posrescale

    beam = get_gauss_beam(fwhm,pixsize,band)
    sim_map = convolve(cmap, beam, boundary='wrap')

    if superplot:
        plt.imshow(outmap,extent=(0,300,0,300),clim=[0.0,0.15],origin=0)
        plt.colorbar()
        plt.title('clus_popmap: Lensed map')
        plt.show()

    if savemaps:
        hdx = fits.PrimaryHDU(maps['signal'],maps['shead'])
        sz = fits.PrimaryHDU(outmap,hdx.header)
        if os.path.isfile(config.SIMBOX + 'lensedmap_' + name + '_' + band + '.fits'):
            os.remove(config.SIMBOX + 'lensedmap_' + name + '_' + band + '.fits')
        sz.writeto(config.SIMBOX + 'lensedmap_' + name + '_' + band + '.fits')

    return sim_map

def get_gauss_beam(fwhm, pixscale, band, nfwhm=5.0, oversamp=1):
    # print('get_gauss_beam: fwhm',fwhm)
    # print('get_gauss_beam: pixscale',pixscale)
    retext = round(fwhm * nfwhm / pixscale)
    if retext % 2 == 0:
        retext += 1

    bmsigma = fwhm / math.sqrt(8 * math.log(2))

    beam = Gaussian2DKernel(bmsigma / pixscale, x_size=retext,
                            y_size=retext, mode='oversample',
                            factor=oversamp)
    beam *= 1.0 / beam.array.max()
    # plt.imshow(beam)
    # plt.savefig('beam_%s.png' %(band))
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
    pixsize = [6.0, 8.33333, 12.0]
    maps, err = clus_get_data(clusname=clusname, resolution=resolution, bolocam=bolocam, verbose=verbose)
    band = maps[0]['band']
    ltfile = config.SIMBOX + 'a0370' + '_image_' + band + '.dat'
    # ltfile = '/data/zemcov/clusters/sandbox/a0370_image_PSW.dat'
    clus_popmap(ltfile,maps[0],band,'a0370',pixsize[0],superplot=1,savemaps=1)
