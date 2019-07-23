
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
from clus_get_lambdas import *
from clus_get_relsz import *


def add_sziso(maps,yin,tin,
              verbose = 0):
    errmsg = False
#   Now the bolocam data is in the SPIRE format
#   we can do a loop over the map adding in the false sz signal
    #This code currently doesnt work

    mapsize = len(maps)

    if verbose:
        print('Fetching BOLOCAM maps')

    data_file = str('bolocam/data/' + maps[0]['name'] + '.sav')
    data_dict = scipy.io.readsav(CLUSDATA + data_file,python_dict = True)
    cluster_struct = list(data_dict.values())
    # Since the bolocam data is in a slightly differnt format we need to convert it to
    # the same format as the hermes data
    bolocam,err = clus_convert_bolocam(cluster_struct,verbose=verbose)

    # See how to incorporate this error message
    # if not ccb_sucess:
    #     errmsg = str('clus_convert_bolocam exited with error: ' + ccb_errmsg)
    #     if verbose:
    #         return None, errmsg
            # Has a go to the err handerler right here in idl

#   Set this to the size of bolocam.deconvolved image
    sziso = fits.PrimaryHDU()
    naxis = bolocam[0]['deconvolved_image'][0].shape
    # print('naxis =',naxis)
    # exit()
#   MKHDR replacement
    # sziso.header = bolocam[0]['deconvolved_image']
    temphead = sziso.header
#   Need to ask what this does
    crpix1 = int(naxis[0] / 2)
    crpix2 = int(naxis[1] / 2)

#   this needs to be put into a temptemphead for the bolocam stuff
#   What is the mapts eleent thats getting used instead
#   use astropy to append to FITS temphead, i.e. hdr.append('CRVAL1')
    # temphead.append('CRVAl1', 'CRVAL2', 'CRPIX1', 'CRPIX2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2',
    #                 'EPOCH', 'EQUINOX', 'CTYPE1', 'CTYPE2')
    print('bolo data',-bolocam[0]['deconvolved_image_resolution_arcmin'][0])
    temphead.set('CRVAL1' , bolocam[0]['deconvolved_image_ra_j2000_deg'][0][crpix1,crpix2])
    temphead.set('CRVAL2' , bolocam[0]['deconvolved_image_dec_j2000_deg'][0][crpix1,crpix2])
    temphead.set('CRPIX1' , crpix1)
    temphead.set('CRPIX2' , crpix2)
    # Why is this one negative???
    temphead.set('CD1_1' , -bolocam[0]['deconvolved_image_resolution_arcmin'][0] / 60.0)
    temphead.set('CD1_2' , 0)
    temphead.set('CD2_1' ,  0)
    temphead.set('CD2_2' , bolocam[0]['deconvolved_image_resolution_arcmin'][0] / 60.0)
    temphead.set('EPOCH' , 2000.00)
    temphead.set('EQUINOX' , 2000.00)
    temphead.set('CTYPE1' , 'RA---TAN')
    temphead.set('CTYPE2' , 'DEC---TAN')

    x = np.arange(naxis[0]) - 14.5
    y = np.arange(naxis[1]) - 14.5
    naxis = np.array(naxis)
    rad = np.zeros((naxis[0],naxis[1]))

    for ix in range(naxis[1]):
        rad[:,ix] = (np.tile(x[ix],naxis[0])**2+y**2)**(1/2)
    '''
    Add in section to find all values in rad that are greater than 10
    Loop through nested for loops because rad is multi dimensional
    Called "outer" , use in szmap
    '''
    # Outer marks the indices of the outside of the circle used to add in the sz effect
    n = 0
    outer = []
    for i in range(len(rad)):
        for j in range(len(rad)):
            if rad[i,j] > 10:
                outer.append(n)
            n +=1
#   right syntax???
    szmap = -1 * bolocam[0]['deconvolved_image'][0]
    szmap = (np.array(szmap)).flatten()
    szmap = szmap - np.mean(szmap[outer])
    szmap = szmap / max(szmap)
    '''Up until this section is working correctly'''
    for imap in range(mapsize):
# Below are the changes to units of Jy to beam
# calculate the surface brightness to put in each pixel at this sim step
# Need to fix maps right now it is written as a pointer
        if imap == 1:
            nu = 3e5 / CLUS_GET_LAMBDAS((maps[imap]).band)
            print(nu)
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
    exit()
    return maps


if __name__ == '__main__':
    add_sziso(norm=1,factor = 0)
