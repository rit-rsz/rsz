
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
import sys
sys.path.append('../utilities')
from clus_convert_bolocam import *
from clus_get_lambdas import *
from config import *
sys.path.append('../sz')
from clus_get_relsz_v2 import *
from astropy.wcs import WCS
from astropy import units as u



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

    for imap in range(mapsize):
# Below are the changes to units of Jy to beam
# calculate the surface brightness to put in each pixel at this sim step
# Need to fix maps right now it is written as a pointer
        # Only on the PLW Band
        if imap == 2:
            # This is the work around until the config is fixed
            yin_coeff = [2.50,1.91,2.26,3.99,1.36,2.42,1.59,1.90,3.99]
            yin = [x*1e-4 for x in yin_coeff]

            tin = [7.2,10.1,7.7,9.8,4.5,8.6,7.8,5.5,10.9]
            nu = 3e5 / clus_get_lambdas((maps[imap]['band']))
            # Change name from sz_wrapper to name of file
            dI_raw,errmsg = sz_wrapper(nu,y=yin) #te=tin
            if errmsg:
                if verbose:
                    new_errmsg = 'Clus_get_relSZ exited with error'+errmsg
                return None, new_errmsg
            dI = dI_raw * 1e26
    #       print,CLUS_GET_LAMBDAS((*maps[imap]).band),dI
    #       this accounts for a conversion to surface brightness in a later step
    #       dI = dI / (1.13d0 * (60./3600.)^2 * (!DTOR)^2 * 1e9)
            #np.tile instead of idls replicate
            # print(maps[imap]['widtha'])
            # print(dI)
            '''Up until this section is working correctly'''
            '''Can Probably still use tile, after the test with zeros it was show the error wa with
                incorporating the mulitplication to szmap. szmap is I think a much larger size and
                it is having a hard time combineig the two'''
            dI_converted = dI*((1.13*(maps[imap]['widtha'])/3600.)**2 * (pi/180)**2)
            # szin = (np.tile(x,(naxis[0],naxis[1]))) * szmap
            szin = dI_converted * szmap
            szin = np.reshape(szin,(naxis[0],naxis[1]))

    #       Have to interpolate to SPIRE map size
            # using hcongrid from astropy to replace HASTROM
            '''Need to writefits szin inorder to use this function'''
            hdu = fits.PrimaryHDU(szin,temphead)
            # hdu.header = temphead
            print(hdu.header['CRVAL1'])
            exit()
            # hdul = fits.HDUList([hdu])
            # Use for debugging
            # hdul.writeto('new1.fits')
            w = WCS(hdu.header)
            # ra = np.array(cats[i]['sra']) * u.deg
            # dec = np.array(hdu) * u.deg
            c = SkyCoord(hdu.header['CRVAL1'],hdu.header['CRVAL2'])
            px, py = skycoord_to_pixel(c, w, 1)
            wc = w.all_world2pix(temphead['CRVAL1'],temphead['CRVAL2'])
            print(wc)
            exit()
            szin.writeto('add_sz_%s.fits' % (maps[imap]['name']))
            w = WCS('add_sz_%s.fits' % (maps[imap]['name']))
            szinp = w.world2pix(naxis[0],naxis[1])
            print(szinp)
            exit()
            maps[imap]['signal'] = maps[imap]['signal'] + szinp



#       Need to then set something for maps thats back to the dictonary format
    exit()
    return maps


if __name__ == '__main__':
    add_sziso(norm=1,factor = 0)
