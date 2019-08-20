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
# matplotlib.use('TkAgg')

def compute_rings(maps, params, binwidth, superplot=1, verbose=1, noconfusion=None):

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
        if verbose:
            print('Average for ' + maps[m]['band'])

        # make sizes for different radial averages
        mapsize = maps[m]['astr']['NAXIS']
        pixsize = maps[m]['pixsize']
        maxrad = ceil(pixsize * np.amax(mapsize) / sqrt(2.0))
        print('Mapsize: ',mapsize,'Pixsize: ',pixsize,'Maxrad: ',maxrad,'nbins: ',nbins)
        # make array objects to fill later
        midbinp = np.empty(nbins)
        midwei = np.empty(nbins)
        fluxbin = np.empty(nbins)
        hitbin = np.empty(nbins)
        errbin = np.empty(nbins)

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
        px, py = skycoord_to_pixel(c, w, origin=1)
        # ===============================================================

        # retrieve noise mask and set base noise level of map
        conv = pixsize / binwidth
        ''' Not sure if this is the right calfac'''
        calfac = (1*10**-6) / (1.13*25.0**2*( pi /180.0)**2/3600.0**2)
        confnoise = confusionnoise[m]
        mask = make_noise_mask(maps, m)
        # ===============================================================
        ''' This is currently where i'm working, this whole section is f*d up'''
        # not exactly sure what this is doing but looks like radius of bin rings
        tempmap = np.empty((int(mapsize[0]), int(mapsize[1])))

        '''more verbose but easier to debug'''
        ''' yeah I could... just seemed a bit verbose'''
        for i in range(mapsize[0]):
            for j in range(mapsize[1]):
                # nthisrad = 0 # keep a counter for how many times we find a minimum
                maps[m]['mask'][i,j] = maps[m]['mask'][i,j] + mask[i,j]
                thisrad = pixsize * sqrt((i - px)**2+(j - py)**2)
                tempmap[i,j] = thisrad
                midbin_fill = [abs(thisrad - x) for x in midbin]
                for k in range(len(midbin)-1):
                    if midbin_fill[k] == np.min(midbin_fill):
                        # print('#################################')
                        # print(k)
                        midbinp[k] = midbinp[k] + thisrad
                        midwei[k] = midwei[k] + 1
                        if  maps[m]['mask'][i,j] == 0 and (k <= maxrad) :
                            thiserr = calfac * sqrt(maps[m]['error'][i,j]**2 + confnoise**2)
                            # print('THISERR: ',thiserr,'CONFNOISE: ',confnoise,'CALFAC: ',calfac)
                            # fluxbin[k] = (calfac * maps[m]['srcrm'][i,j] / thiserr**2)
                            fluxbin[k] = fluxbin[k] + (calfac / thiserr**2)
                            hitbin[k] = hitbin[k] + 1.0 / thiserr**2

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
            if hitbin[i] > 0 :
                fluxbin[i] = fluxbin[i] / hitbin[i]
                errbin[i] = sqrt(1.0 / hitbin[i])
            else :
                midbinp[i] = np.nan
                fluxbin[i] = np.nan
                errbin[i] = np.nan

        print('FLUXBIN: ',fluxbin)
        print('MIDBIN: ',midbinp)
        exit()
        #save new bin data to dictionary & return to fitsz
        radave = [None]*ncols
        radave[m] = {'band' : maps[m]['band'],
                     'midbin' : midbinp,
                     'fluxbin' : fluxbin,
                     'errbin' : errbin}

        # print(radave[m]['midbin'])

        if superplot:
            plt.scatter(midbinp,fluxbin)
            plt.title('Clus Compute Rings: Radial Averages for %s' %(maps[m]['band']))
            plt.xlabel('Radius (arcsec)')
            plt.ylabel('Signal (MJy/sr)')
            plt.xlim((0,600))
            plt.ylim((-1000,1000))
            plt.show()


    return radave

if __name__ == '__main__' :
    maps,err = get_data('a0370')
    params,err = get_clus_params('a0370')
    compute_rings(maps,params,30.0)
