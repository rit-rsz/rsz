
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
from math import *
from astropy.io import fits
import os
# import sys
# sys.path.append('../')
from config import *
import sys
sys.path.append('../utilities')
from clus_convert_bolocam import *


def add_sziso(maps,yin,tin,
              verbose = 0):
    errmsg = False
#   Now the bolocam data is in the SPIRE format
#   we can do a loop over the map adding in the false sz signal
    #This currently doesnt work
    mapsize = len(maps)

    if verbose:
        print('Fetching BOLOCAM maps')

#   Need to check how its calling the maps in simmaps

    # This is needed until get_simmaps is fixed to out put maps[names] as a str
    clusname = maps[0]['name'].astype(str)
    # print(a[0])
    data_file = str('bolocam/data/' + clusname[0] + '.sav')
    print(data_file)
    data_dict = scipy.io.readsav(CLUSDATA + data_file,python_dict = True)
    # example of how you need to call this list
    # data_dict.setflags(write=1)
    cluster_struct = list(data_dict.values())
    print(cluster_struct['deconvolved_image'])

    exit()
    bolocam = clus_convert_bolocam(cluster_struct,verbose=verbose)
    exit()
    if not ccb_sucess:
        errmsg = str('clus_convert_bolocam exited with error: ' + ccb_errmsg)
        if verbose:
            return None, errmsg
            # Has a go to the err handerler right here in idl

#   Set this to the size of bolocam.deconvolved image
    # naxis =
#   Some kind of temphead maker??? Fuction??
#   MKHDR
    temphead = fits.Header(bolocam['deconvolved_image'])
#   Need to ask what this does
    crpix = naxis[1] / 2
    crpix = naxis[2] / 2
#   this needs to be put into a temptemphead for the bolocam stuff
#   What is the mapts eleent thats getting used instead
#   use astropy to append to FITS temphead, i.e. hdr.append('CRVAL1')
    temphead.append('CRVAl1', 'CRVAL2', 'CRPIX1', 'CRPIX2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2',
                    'EPOCH', 'EQUINOX', 'CTYPE1', 'CTYPE2')
    temphead['CRVAL1'] = [bolocam['deconvolved_image_ra_j2000_deg']['crpix1'],bolocam['deconvolved_image_ra_j2000_deg']['crpix2']]
    temphead['CRVAL2'] = [bolocam['deconvolved_image_dec_j2000_deg']['crpix1'],bolocam['deconvolved_image_dec_j2000_deg']['crpix2']]
    temphead['CRPIX1'] = crpix1
    temphead['CRPIX2'] = crpix2
    # Why is this one negative???
    temphead['CD1_1'] = -1 * bolocam['deconvolved_image_resolution_arcmin'] / 60.0
    temphead['CD1_2'] = 0
    temphead['CD2_1'] = 0
    temphead['CD2_2'] = bolocam['deconvolved_image_resolution_arcmin'] / 60.0
    temphead['EPOCH'] = 2000.00
    temphead['EQUINOX'] = 2000.00
    temphead['CTYPE1'] = 'RA---TAN'
    temphead['CTYPE2'] = 'DEC---TAN'

#   This was findgen in idl. I may have to chage this to np.zeros I need to test the output in idl to confirm
    x = np.arange(naxis[1]) - 14.5
    y = np.arange(naxis[2]) - 14.5

    rad = np.empty(naxis[1],naxis[2])


#   Need to figure out how to use replicate in a fashion for python
    for ix in range(naxis[2]):
        rad[:,ix] = sqrt(replicate())

    circ = where(rad < 10)

#   right syntax???
    szmap = -1 * bolocam.deconvolved_image
    szmap = szmap - mean(szmap[outer])
    szmap = szmap / max(szmap)
    print('hello')
    for imap in range(mapsize):
# Below are the changes to units of Jy to beam
# calculate the surface brightness to put in each pixel at this sim step
# Need to fix maps right now it is written as a pointer
        nu = 3e5 / CLUS_GET_LAMBDAS((maps[imap]).band)
        dI = 1e26*clus_get_relSZ(nu,yin,tin)
#       print,CLUS_GET_LAMBDAS((*maps[imap]).band),dI
#       this accounts for a conversion to surface brightness in a later step
#       dI = dI / (1.13d0 * (60./3600.)^2 * (!DTOR)^2 * 1e9)
        #np.tile instead of idls replicate
        szin = repeat(dI*((1.13*(maps[imap]).widtha/3600.)^2 * \
                    (pi/180)^2),naxis[1],naxis[2]) * szmap

#       Have to interpolate to SPIRE map size
        HASTROM

#       Need to then set something for maps thats back to the dictonary format

    return maps


if __name__ == '__main__':
    add_sziso(norm=1,factor = 0)
