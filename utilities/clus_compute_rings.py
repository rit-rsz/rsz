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
from clus_make_noise_mask import *
from clus_get_data import clus_get_data
from clus_get_clusparams import *
from astropy.wcs.utils import skycoord_to_pixel
import config
# from FITS_tools.load_header import load_header
import matplotlib.pyplot as plt
import matplotlib
# matplotlib.use('TkAgg')

def clus_compute_rings(maps, params, binwidth, superplot=0, verbose=1, noconfusion=None):
    # init params
    bands = ['PSW', 'PMW', 'PLW']
    ncols = len(maps)
    radave = [None]*ncols

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
    for m in range(ncols):
        clusname = maps[m]['name']
        srcrm = maps[m]['signal']
        print(maps[m]['band'])
        plt.imshow(srcrm)
        plt.show()
        if verbose:
            print('Average for ' + maps[m]['band'])

        # make sizes for different radial averages
        mapsize = maps[m]['astr']['NAXIS']
        pixsize = maps[m]['pixsize']
        maxrad = ceil(pixsize * np.amax(mapsize) / sqrt(2.0))
        if verbose:
            print('Mapsize: ',mapsize,'Pixsize: ',pixsize,'Maxrad: ',maxrad,'nbins: ',nbins)
        # make array objects to fill later
        midbinp = np.zeros(nbins)
        midwei = np.zeros(nbins)
        fluxbin = np.zeros(nbins)
        hitbin = np.zeros(nbins)
        errbin = np.zeros(nbins)

        # make radbin
        step = maxrad / (nbins-1) + 3.0 * float(m) #double check this to convince our selves this is what we want.
        radbin = np.arange(0.0,maxrad+step,step)

        # make midbin
        bin_size = np.absolute([x / 2.0 for x in np.diff(radbin,1)])
        midbin = np.append([x + bin_size[0] for x in radbin],maxrad)

        #convert RA/DEC to pixel coordinates
        ra = params['fidrad'] * u.deg
        dec = params['fidded'] * u.deg
        c = SkyCoord(ra, dec)
        # grabbing header from maps file
        hdul = fits.open(maps[m]['file']) #ideally we would be able to use shead here instead.
        w = wcs.WCS(hdul[1].header)
        #converting ra/dec to pixel coords
        px, py = skycoord_to_pixel(c, w, origin=0)
        # ===============================================================

        # retrieve noise mask and set base noise level of map
        conv = pixsize / binwidth
        confnoise = confusionnoise[m]
        mask = clus_make_noise_mask(maps, m)
        # ===============================================================
        # once filled, shows the radial bins from thisrad
        tempmap = np.zeros((int(mapsize[0]), int(mapsize[1])))

        # find the flux that falls closest a given radius (thisrad) for a given pixel (ipix,jpix)
        for ipix in range(0, mapsize[0]-1):
            for jpix in range(0, mapsize[1]-1):
                maps[m]['mask'][ipix,jpix] = maps[m]['mask'][ipix,jpix] + mask[ipix,jpix]
                thisrad = pixsize * sqrt((ipix - py)**2+(jpix - px)**2)
                tempmap[ipix,jpix] = thisrad

                midbin_fill = [abs(thisrad - x) for x in midbin]
                rad = np.where(midbin_fill == np.min(midbin_fill))
                rad = rad[0][0]
                if rad < nbins :
                    midbinp[rad] = midbinp[rad] + thisrad
                    midwei[rad] = midwei[rad] + 1
                    if maps[m]['mask'][ipix,jpix] == 0 and (rad <= maxrad) :
                        #calculating our value for sigma^2
                        thiserr = maps[m]['calfac'] * sqrt(maps[m]['error'][ipix,jpix]**2 + confnoise**2)
                        #summing up the flux * 1 / sigma^2
                        '''
                            srcrm is the new pcat_spire map, which may already be in units of
                            mJy/pixel rather than Jy/beam like the other maps
                        '''
                        fluxbin[k] = fluxbin[k] + (maps[m]['calfac'] * maps[m]['srcrm'][i,j] / thiserr**2)
                        #summing up 1 / sigma^2
                        hitbin[rad] = hitbin[rad] + 1.0 / thiserr**2
                        ''' for testing without real subtracted map
                        thiserr = 1
                        fluxbin[rad] = fluxbin[rad] + ( maps[m]['signal'][ipix, jpix] / thiserr**2)
                        '''

        # =========================================================================================
        # file = config.FITSOUT + 'radmap_' + bands[m] + '_' + maps[m]['name'] + '.fits'
        # if os.path.isfile(file):
        #     os.remove(file)
        # x = load_header(str(maps[m]['shead']))
        # hdu = fits.PrimaryHDU(tempmap,x)
        # hdu.writeto(file)
        # ===========================================================================================

        #finding averages and changing some value to nan if error points to that.
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

            fig = plt.figure()
            new_midbin = [(x/pixsize) + px for x in midbinp]

            for k in range(nbins):
                circle = plt.Circle((px,py),new_midbin[k]-px,fill=False)
                plt.gca().add_artist(circle)
            plt.imshow(tempmap, alpha=0.9)
            plt.scatter(new_midbin,[py]*nbins,c='r')
            plt.title('Midbins for %s' %(maps[m]['band']))
            plt.show()

    return radave


if __name__ == '__main__' :
    maps,err = clus_get_data('a0370')
    params,err = clus_get_clusparams('a0370')
    clus_compute_rings(maps,params,30.0)
