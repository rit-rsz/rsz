################################################################################
# NAME : get_simmaps.py
# DATE STARTED : June 18, 2019
# AUTHORS : Benjamin Vaughan & Dale Mercado
# PURPOSE : Fetches the simmaps
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS : clusname - name of the cluster
#         nsim - name of the sim number
#         simflag - 0 not using sims
#                   1 100 generation sims CLuster_simulation
#                   2 200 generation sims Benthermin_simulations
#         sb        Flag for surface Brightness
#         xc        Flag for Cross Cleaned
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
        # print(map.shead[0])
        # print(map.shead[0].CD1_1)
        head = fits.Header()
        map.shead[0] = map.shead[0].astype(str)
        for i in range(map.shead[0].shape[0]):
            line = map.shead[0][i]
            split_line = line.split('=')
            key = split_line[0]
            if len(split_line) == 2:
                split_line = split_line[1].split('/')
                value = split_line[0]
                comment = split_line[1]
                head.set(key, value, comment)
        map.ehead[0] = map.ehead[0].astype(str)
        ehead = fits.Header()
        for i in range(map.ehead[0].shape[0]):
            line = map.ehead[0][i]
            split_line = line.split('=')
            key = split_line[0]
            if len(split_line) == 2:
                split_line = split_line[1].split('/')
                value = split_line[0]
                comment = split_line[1]
                ehead.set(key, value, comment)
        astr = fits.Header()
        # print(map.astr[0])
        # map.astr[0] = map.astr[0].astype(str)
        # for i in range(map.astr[0].shape[0]):
        #     line = map.astr[0][i]
        #     split_line = line.split('=')
        #     print(split_line)
        #     key = split_line[0]
        #     if len(split_line) == 2:
        #         split_line = split_line[1].split('/')
        #         value = split_line[0]
        #         comment = split_line[1]
        #         astr.set(key, value, comment)

        new_map = {'name' : map['name'].astype(str)[0],
                   'file' : map['file'][0].astype(str),
                   'band' : map.band.astype(str)[0],
                   'signal' : map.signal[0].astype(float),
                   'srcrm' : map.srcrm[0],
                   'xclean' : map.xclean[0],
                   'error' : map.error[0],
                   'mask' : map.mask[0],
                   'shead' : head,
                   'ehead' : ehead,
                   'astr' : map.astr[0],
                   'pixsize' : map.pixsize[0].astype(float),
                   'psf' : map.psf[0].astype(float),
                   'width' : map.width[0].astype(float),
                   'widtha' : map.widtha[0].astype(float),
                   'calfac' : map.calfac[0].astype(float),
                   'jy2mjy' : map.jy2mjy[0].astype(float)}
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
