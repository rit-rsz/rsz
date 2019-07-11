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


def get_simmaps(clusname, nsim, simflag=1, sb=0, xc=0, verbose=0, maps = None):

    band = ['PSW', 'PMW', 'PLW'] #setting the different bands
    bigpix = [6.0,8.3,12.0] #arcsec/pixel for each of the bands.
    nbands = len(band)

    if simflag == 0 :
        return maps

    elif simflag == 1: #if we want to get simulated data.
        maps = [] #initializing maps.
        for i in range(nbands):
            mapfilename = config.CLUSDATA + 'sim_clusters/' + clusname + '_' + band[i] + '.sav'
            thismap = sp.readsav(mapfilename,python_dict = True) #extracting a python dictionary to contain the data for each map from a .sav file
            #have to use format thismap['thismap'] because of the way data gets extracted from .sav file in python.
            #also thismap['thismap'] is a numpy record array.
            thismap['thismap'].calfac = 1.0 / (thismap['thismap'].calfac * (thismap['thismap'].pixsize)* thismap['thismap'].pixsize)
            if sb:
                overmap = fits.open(config.CLUSDATA + 'sim_clusters/' + clusname +
                '_' + band[i] + '_sb.fits') #opening a fits file to get some data.
                overmap = overmap[0].data
                thismap['thismap'].xclean[0] = overmap
                #this block of code is commented out because the assignment destination is read-only.
                # whpl = []
                # for j in range(overmap.shape[0]):
                #     for k in range(overmap.shape[1]):
                #         if np.isfinite(overmap[j,k]) == False:
                #             whpl.append([j,k])
                # thismap['thismap'].mask[0][whpl] = 1
            elif xc:
                overmap = fits.open(config.CLUSDATA + 'sim_clusters/' + clusname +
                '_' + band[i] + '_xc.fits') #opening a fit file to get some data.
                overmap = overmap[0].data
                thismap['thismap'].xclean[0] = overmap
                #this block of code is commented out because the assignment destination is read-only.
                # for j in range(overmap.shape[0]):
                #     for k in range(overmap.shape[1]):
                #         if np.isfinite(overmap[j,k]) == False:
                #             whpl.append([j,k])
                # thismap['thismap'].mask[0][whpl] = 1
            maps.append(thismap['thismap'])
    else:
        maps = [] #initializing maps.
        for i in range(nbands):
            # thissim = lead_zero(nsim)
            thissim = str(nsim)
            mapfilename = config.CLUSDATA + 'bethermin_sims/' + clusname +'/' + clusname + '_' + band[i] + '_sim0' + thissim + '.sav'
            thismap = sp.readsav(mapfilename,python_dict = True) #get data from a .sav file.

            #creating noise maps and signal maps.
            mapsize = thismap['thismap'].signal[0].shape
            noisemap = 1.0 * thismap['thismap'].error[0] * np.random.standard_normal((mapsize[0], mapsize[1]))
            thismap['thismap'].signal[0] = thismap['thismap'].signal[0] + noisemap
            maps.append(thismap['thismap'])

    new_maps = []
    for map in maps:
        # print(map.dtype.names)
        new_map = {'name' : map['name'].astype(str)[0],
                   'file' : map['file'][0].astype(str),
                   'band' : map.band.astype(str)[0],
                   'signal' : map.signal[0],
                   'srcrm' : map.srcrm[0],
                   'xclean' : map.xclean[0],
                   'error' : map.error[0],
                   'mask' : map.mask[0],
                   'shead' : map.shead[0],
                   'ehead' : map.ehead[0],
                   'astr' : map.astr[0],
                   'pixsize' : map.pixsize[0],
                   'psf' : map.psf[0],
                   'width' : map.width[0],
                   'withda' : map.widtha[0],
                   'calfac' : map.calfac[0],
                   'jy2mjy' : map.jy2mjy[0]}
        new_maps.append(new_map)
    return new_maps ,None

# def lead_zero(nsim):
#     nsim = str(nsim)
#     if len(nsim) == 3 :
#         out = "%01d" %int(nsim)
#     elif len(nsim) == 2 :
#         out = "%02d" %int(nsim)
#     elif len(nsim) == 1 :
#         # out = f"{1:03d}"
#         out = "%03d" %int(nsim)
#     return out

if __name__ == "__main__":
    maps = get_simmaps('a0370', nsim=200, simflag=2)
