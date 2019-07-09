################################################################################
# NAME : get_simmaps.py
# DATE STARTED : June 18, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE : Fetches the simmaps
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS : clusname - name of the cluster
#         nsim - name of the sim number
#         simflag - 1 = do something
#                   0 = do something else
#         sb      - 1 = do something
#                 - 0 = do something else
#         xc      - 1 = do something
#                 - 0 = do something else
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import numpy as np
from astropy.io import fits
import random
import scipy.io as sp
import config


def get_simmaps(clusname, nsim, simflag=1, sb=0, xc=0, verbose=0):

    if simflag > 0 and nsim == 0:
        print('simmap set but nsim not supplied! Aborting.')
        exit()
    cols = ['PSW', 'PMW', 'PLW']
    bigpix = [6.0,8.3,12.0] #arcsec/pixel
    ncols = len(cols)
    maps = []
    for i in range(ncols):
        maps.append([])
    if simflag == 1:
        for i in range(ncols):
            mapfilename = config.CLUSDATA + 'sim_clusters/' + clusname + '_' + cols[i] + '.sav'
            thismap = sp.readsav(mapfilename, python_dict = True)
            print(thismap['thismap'].calfac)
            thismap['thismap'].calfac = 1.0 / (thismap['thismap'].calfac * (thismap['thismap'].pixsize)* thismap['thismap'].pixsize)
            if sb:
                overmap = fits.open(config.CLUSDATA + 'sim_clusters/' + clusname +
                '_' + cols[i] + '.fits') #replaced _sb.fits with .fits because there was no _sb.fits file?
                overmap = overmap[1].data
                thismap['thismap'].xclean[0] = overmap
                whpl = np.where(np.isfinite(overmap) == False)
                for index in whpl:
                    thismap['thismap'].mask[0][index] = 1
            elif xc:
                overmap = fits.open(config.CLUSDATA + 'sim_clusters/' + clusname +
                '_' + cols[i] + '.fits') #replaced _xc.fits with .fits because no _xc.fits file?
                overmap = overmap[1].data
                thismap['thismap'].xclean[0] = overmap
                whpl = np.where(np.isfinite(overmap) == False)
                for index in whpl:
                    thismap['thismap'].mask[0][index] = 1 #read only array ???
                thismap['thismap'].mask = new_arr
            maps[i] = thismap['thismap']
    else:
        for i in range(ncols):
            thissim = str(nsim) #idk what format = I04 does...
            mapfilename = config.CLUSDATA + 'bethermin_sims/' + clusname +'/' + clusname + '_' + cols[i] + '_sim0' + thissim + '.sav'
            thismap = sp.readsav(mapfilename)
            mapsize = thismap['thismap'].signal[0].shape #still don't know what thismap is.
            noisemap = 1.0 * thismap['thismap'].error[0] * np.random.standard_normal((mapsize[0], mapsize[1]))
            thismap['thismap'].signal[0] = thismap['thismap'].signal[0] + noisemap #STILL don't know what thismap is
            maps[i] = thismap['thismap']
    return maps ,None

if __name__ == "__main__":
    maps = get_simmaps('rxj1347', 100, simflag=1, xc=1)
