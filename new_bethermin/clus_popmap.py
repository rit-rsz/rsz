
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
from astropy.convolution import convolve_fft
from PIL import Image
from datetime import datetime
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

    if superplot:
        plt.scatter(ra,dec,s=2,c=flux)
        plt.colorbar()
        plt.title('Start of Clus Popmap')
        plt.show()

    header = maps['shead']
    wcs = WCS(header)
    coords = SkyCoord(ra=ra*u.deg,dec=dec*u.deg)
    x,y = skycoord_to_pixel(coords, wcs)

    if superplot:
        plt.scatter(x,y,s=2,c=flux)
        plt.colorbar()
        plt.title('clus_popmap: after skycoord_to_pixel')
        plt.show()

    x_size = maps['signal'].shape[0]
    y_size = maps['signal'].shape[1]
    # x_size = 300
    # y_size = 300
    outmap = np.zeros((x_size,y_size))
    num = 0
    fig, axs = plt.subplots(1, 2)

    # print(datetime.now(), 'Time before manually multiplying through gaussians')
    for i in range(len(flux)):
        if x[i] > 0 and x[i] < y_size and y[i] > 0 and y[i] < x_size:
            kern = makeGaussian(y_size,x_size, fwhm = fwhm/pixsize, center=(int(x[i]),int(y[i])))
            kern = kern / np.max(kern)
            norm = flux[i]
            psf = kern * norm
            outmap = outmap + psf
        else :
            num += 1
    # axs[0].imshow(outmap)
    #
    # print(datetime.now(), 'Time after manually multiplying through gaussians')

    # emptymap = np.zeros((x_size, y_size))
    # img = Image.fromarray(emptymap)
    # for i in range(len(flux)):
    #     img.putpixel((x[i],y[i]), flux[i])


    # im = np.asarray(img)
    #
    # kern = makeGaussian(15, 15, fwhm= fwhm / pixsize)
    # kern = kern / np.max(kern)
    # kern = padGaussian(im, kern)
    #
    # outmap = convolve_fft(outmap, kern)
    #
    # axs[1].imshow(outmap)

    print('number of sources outside map: ',num)

    if superplot:
        plt.imshow(outmap,extent=(0,300,0,300),clim=[0.0,0.15],origin=0)
        plt.colorbar()
        plt.title('clus_popmap: Lensed map')
        plt.show()

    if savemaps:
        hdx = fits.PrimaryHDU(maps['signal'],maps['shead'])
        sz = fits.PrimaryHDU(outmap,hdx.header)
        sz.writeto(config.SIMBOX + 'lensedmap_' + name + '_' + band + '.fits',overwrite=True)

    return outmap

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
    clus_popmap(ltfile,maps[0],band,'a0370',pixsize[0], 18, superplot=1,savemaps=1)
