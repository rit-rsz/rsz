################################################################################
# NAME : compute_rings.py
# DATE STARTED : June 24, 2019
# AUTHORS : Benjamin Vaughan & Victoria Butler
# PURPOSE : computes the radial average
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
from math import *
import numpy as np
from astropy import wcs
from astropy.io import fits
from astropy.coordinates import SkyCoord
import sys, os
from astropy import units as u
sys.path.append('../utilities')
from make_noise_mask import make_noise_mask
from get_data import get_data
from get_clusparams import get_clus_params
from astropy.wcs.utils import skycoord_to_pixel
import config
# from FITS_tools.load_header import load_header
import matplotlib.pyplot as plt
import matplotlib
import time
# matplotlib.use('TkAgg')

def compute_rings(maps, params, binwidth, superplot=1, verbose=1, noconfusion=None):
    t0 = time.time()
    # init params
    bands = ['PSW', 'PMW', 'PLW']
    ncols = len(maps)
    if noconfusion:
        confusionnoise = np.empty(ncols)
    else:
        ar = [5.8, 6.3, 6.8]
        confusionnoise = np.array([x*0.001 for x in ar])

    # for exclusively calculating largest nbins
    nbins = 0
    for i in range(ncols):
        mapsize = maps[i]['astr']['NAXIS']
        pixsize = maps[i]['pixsize']
        maxrad = ceil(pixsize * np.amax(mapsize) / sqrt(2.0))
        nbinsp = int(maxrad / binwidth) + 1
        if nbinsp > nbins :
            nbins = int(nbinsp)

    # create bins for different metrics (mid,flux,err,rad,hit,max)
    radave = []
    for m in range(ncols):
        clusname = maps[m]['name']
        # hdul = fits.open('/home/butler/bitten/SPIRE/cluster_analysis/plots/clus_rings_test_%s_%s.fits' %(bands[m],clusname))
        # # hdul.writeto('test.fits')
        # srcrm = hdul[0].data
        if verbose:
            print('Average for ' + maps[m]['band'])

        # make sizes for different radial averages
        mapsize = maps[m]['astr']['NAXIS']
        pixsize = maps[m]['pixsize']
        maxrad = ceil(pixsize * np.amax(mapsize) / sqrt(2.0))
        print('Mapsize: ',mapsize,'Pixsize: ',pixsize,'Maxrad: ',maxrad,'nbins: ',nbins)
        # make array objects to fill later
        midbinp = np.zeros(nbins)
        midwei = np.zeros(nbins)
        fluxbin = np.zeros(nbins)
        hitbin = np.zeros(nbins)
        errbin = np.zeros(nbins)

        # make radbin
        step = maxrad / (nbins-1) + 3.0 * float(m)
        radbin = np.arange(0.0,maxrad+step,step)

        # make midbin
        bin_size = np.absolute([x / 2.0 for x in np.diff(radbin,1)])
        midbin = np.append([x + bin_size[0] for x in radbin],maxrad)

        #convert RA/DEC to pixel coordinates
        ra = params['fidrad'] * u.deg
        dec = params['fidded'] * u.deg
        c = SkyCoord(ra, dec)
        # grabbing header from maps file
        hdul = fits.open(maps[m]['file'])
        w = wcs.WCS(hdul[1].header)
        #converting ra/dec to pixel coords
        px, py = skycoord_to_pixel(c, w, origin=0)
        # ===============================================================

        # retrieve noise mask and set base noise level of map
        conv = pixsize / binwidth
        confnoise = confusionnoise[m]
        mask = make_noise_mask(maps, m)
        # ===============================================================
        # once filled, shows the radial bins from thisrad
        tempmap = np.zeros((int(mapsize[0]), int(mapsize[1])))

        # find the flux that falls closest a given radius (thisrad) for a given pixel (ipix,jpix)
        for ipix in range(mapsize[0]):
            for jpix in range(mapsize[1]):
                maps[m]['mask'][ipix,jpix] = maps[m]['mask'][ipix,jpix] + mask[ipix,jpix]
                thisrad = pixsize * sqrt((ipix - px)**2+(jpix - py)**2)
                tempmap[ipix,jpix] = thisrad

                midbin_fill = [abs(thisrad - x) for x in midbin]
                rad = np.where(midbin_fill == np.min(midbin_fill))
                rad = rad[0][0]
                if rad < nbins :
                    midbinp[rad] = midbinp[rad] + thisrad
                    midwei[rad] = midwei[rad] + 1
                    if maps[m]['mask'][ipix,jpix] == 0 and (rad <= maxrad) :
                        thiserr = maps[m]['calfac'] * sqrt(maps[m]['error'][ipix,jpix]**2 + confnoise**2)
                        fluxbin[k] = fluxbin[k] + (maps[m]['calfac'] * maps[m]['srcrm'][i,j] / thiserr**2)
                        # fluxbin[rad] = fluxbin[rad] + (maps[m]['calfac'] * srcrm[ipix,jpix] / thiserr**2)
                        hitbin[rad] = hitbin[rad] + 1.0 / thiserr**2

        # =========================================================================================
        # file = config.FITSOUT + 'radmap_' + bands[m] + '_' + maps[m]['name'] + '.fits'
        # if os.path.isfile(file):
        #     os.remove(file)
        # x = load_header(str(maps[m]['shead']))
        # hdu = fits.PrimaryHDU(tempmap,x)
        # hdu.writeto(file)
        # ===========================================================================================

        for i in range(nbins):
            if midwei[i] > 1.0 :
                midbinp[i] = midbinp[i] / midwei[i]
        for j in range(nbins):
            if hitbin[j] > 0 :
                fluxbin[j] = fluxbin[j] / hitbin[j]
                errbin[j] = sqrt(1.0 / hitbin[j])
            else :
                midbinp[j] = np.nan
                fluxbin[j] = np.nan
                errbin[j] = np.nan

        #save new bin data to dictionary & return to fitsz
        radave = [None]*ncols
        radave[m] = {'band' : maps[m]['band'],
                     'midbin' : midbinp,
                     'fluxbin' : fluxbin,
                     'errbin' : errbin}

        if superplot:
            plt.plot(midbinp,fluxbin)
            plt.title('Clus Compute Rings: Radial Averages for %s' %(maps[m]['band']))
            plt.xlabel('Radius (arcsec)')
            plt.ylabel('Signal (MJy/sr)')
            plt.xlim((0,600))
            plt.show()

    return radave


if __name__ == '__main__' :
    maps,err = get_data('a0370')
    params,err = get_clus_params('a0370')
    compute_rings(maps,params,30.0)
