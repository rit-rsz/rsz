
################################################################################
# NAME : add_sziso.py
# DATE STARTED : June 18, 2019
# AUTHORS : Dale Mercado
# PURPOSE : Adding in the SZ isomap
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import scipy.io
import numpy as np
# from config import * #(this line will give me access to all directory variables)
import matplotlib.pyplot as plt
import math
from astropy.io import fits
import os
import pyfits

def add_sziso(maps,YIN=yin,TIN=tin,\
                          verbose = 0,ERRMSG=errmsg,success = 0b):

#   Now the bolocam data is in the SPIRE format
#   we can do a loop over the map adding in the false sz signal
    mapsize = len(maps)

    if verbose:
        print('Fetching Bolocam maps')
#   Need to check how its calling the maps name
    data_file = CLUSDATA + 'bolocam/data/' + maps[0].name + '.sav'
    data = scipy.io.readsav(CLUSDATA + data_file,python_dict = True)
    arr_data = data.values()[0][k][0] # only works on python 2, returns increment
    ''' for python 3 use arr_data = list(data.values())[0][k][0] '''

    bolocam = convert bolocam(cluster_struct,VERBOSE=verbose,\
                                 ERRMSG=ccb_errmsg,SUCCESS=ccb_success)
    if not ccb_sucess:
        errmsg = str('clus_convert_bolocam exited with error: ' + ccb_errmsg)
        if verbose:
            print(errmsg)
            # Has a go to the err handerler right here in idl

#   Set this to the size of bolocam.deconvolved image
    # naxis =
#   Some kind of header maker??? Fuction??
#   MKHDR

#   Need to ask what this does
    crpix = naxis[1] / 2
    crpix = naxis[2] / 2

#   What is the mapts eleent thats getting used instead
    map[0].header['CRVAL1'] = bolocam.deconvolved_image_ra_j2000_deg[crpix1,crpix2]
    map[0].header['CRVAL2'] = bolocam.deconvolved_image_dec_j2000_deg[crpix1,crpix2]
    map[0].header['CRPIX1'] = crpix1
    map[0].header['CRPIX2'] = crpix2
    # Why is this one negative???
    map[0].header['CD1_1'] = -bolocam.deconvolved_image_resolution_arcmin/60
    map[0].header['CD1_2'] = 0
    map[0].header['CD2_1'] = 0
    map[0].header['CD2_2'] = bolocam.deconvolved_image_resolution_arcmin/60
    map[0].header['EPOCH'] = 2000.00
    map[0].header['EQUINOX'] = 2000.00
    map[0].header['CTYPE1'] = 'RA---TAN'
    map[0].header['CTYPE2'] = 'DEC---TAN'

#   This was findgen in idl. I may have to chage this to np.zeros I need to test the output in idl to confirm
    x = np.arange(naxis[1]) - 14.5
    y = np.arange(naxis[2]) - 14.5

    rad = np.zeros(naxis[1],naxis[2])


#   Need to figure out how to use replicate in a fashion for python
    for ix in range(naxis[2]):
        rad[:,ix] = sqrt(replicate())

    circ = where(rad < 10)

#   right syntax???
    szmap = -1 * bolocam.deconvolved_image
    szmap = szmap - mean(szmap[outer])
    szmap = szmap / max(szmap)

    for imap in range(mapsize):
# Below are the changes to units of Jy to beam
# calculate the surface brightness to put in each pixel at this sim step
     nu = 3e5 / CLUS_GET_LAMBDAS((*maps[imap]).band)
     dI = 1e26*clus_get_relSZ(nu,yin,tin)
#    print,CLUS_GET_LAMBDAS((*maps[imap]).band),dI
#     this accounts for a conversion to surface brightness in a later step
#     dI = dI / (1.13d0 * (60./3600.)^2 * (!DTOR)^2 * 1e9)

    szin = replicate(dI*((1.13*(*maps[imap]).widtha/3600.)^2 * $
                    (!DTOR)^2),naxis[1],naxis[2])* szmap)

    
