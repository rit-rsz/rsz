
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
import sys
sys.path.append('../utilities')
from gaussian import makeGaussian
import matplotlib.pyplot as plt
from astropy.io import fits
from FITS_tools.hcongrid import hcongrid
import numpy as np

def clus_popmap(ltfile,maps,band,name,pixsize,loz=None):

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

    naxis = maps['signal'].shape
    crpix1 = int(naxis[0] / 2)
    crpix2 = int(naxis[1] / 2)
    refx = crpix1
    refy = crpix2

    # convert ra/dec to degrees
    print(ra[0:10])
    print(dec[0:10])
    ra = [((-x / 3600.0) + refx) for x in ra]
    dec = [((y / 3600.0) + refy) for y in dec]
    print(ra[0:10])
    print(dec[0:10])

    flux = [10.0**(-k/2.5) for k in mag]
    # if len(loz) == 8 :
    #     flux = [flux,loz['f']]

    # header = maps['shead']
    # wcs = WCS(header)
    # coords = SkyCoord(ra=ra*u.deg,dec=dec*u.deg)
    # x,y = skycoord_to_pixel(coords, wcs, origin=0)
    x = [i / pixsize for i in ra]
    y = [j / pixsize for j in dec]
    print(x[0:10])
    print(y[0:10])

    x_size = maps['signal'].shape[0]
    y_size = maps['signal'].shape[1]
    outmap = np.zeros((x_size,y_size))
    for i in range(len(flux)):
        kern = makeGaussian(y_size,x_size, fwhm = 3, center=(x[i],y[i]))
        kern = kern / np.max(kern)
        norm = flux[i]
        psf = kern * norm
        outmap = outmap + psf
        # plt.imshow(outmap)
        # plt.show()

    # plt.imshow(outmap)
    plt.imsave('/home/butler/rsz/popmap.png',outmap,format='png')
    plt.show()

    return

    # hdx = fits.PrimaryHDU(maps['signal'],maps['shead'])
    # sz = fits.PrimaryHDU(outmap,hdx.header)
    # sz.writeto(config.SIMBOX + 'lensedmap_' + name + '_' + band + '.fits')

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
    band = maps[2]['band']
    ltfile = config.SIMBOX + 'a0370' + '_image_' + band + '.dat'
    clus_popmap(ltfile,maps[2],band,'a0370',pixsize[2])
