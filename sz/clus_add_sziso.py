
################################################################################
# NAME : add_sziso.py
# DATE STARTED : June 18, 2019
# AUTHORS : Dale Mercado
# PURPOSE : Adding in the SZ isomap
# EXPLANATION : Uses bolocam data to create the spectral shape of the sz effect
#               which can then be combined with the peak intensty of the SZ.
#               This is then added into the simmulated 500um map.
# CALLING SEQUENCE :
# INPUTS : maps: Simmulated map objects
#          yin: Compton y paramter values
#          tin: Temperature of electrons
#          verbose: Print out messages
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import scipy.io
import numpy as np
from math import *
import os
import sys
sys.path.append('../utilities')
from clus_convert_bolocam import *
from clus_get_lambdas import *
from config import *
sys.path.append('../sz')
sys.path.append('../source_handling')
from clus_get_data import *
from clus_get_relsz import *
from astropy.io import fits
# from astropy import units as u
from FITS_tools.hcongrid import hcongrid
from IB_model import *  # for testing



def clus_add_sziso(maps,yin,tin,params,ib_test,
              verbose = 0,):
    errmsg = False

    # Now the bolocam data is in the SPIRE format
    # we can do a loop over the map adding in the false sz signal

    mapsize = len(maps)

    if verbose:
        print('Fetching BOLOCAM maps')

    data_file = str('bolocam/data/' + maps[0]['name'] + '.sav')
    data_dict = scipy.io.readsav(CLUSDATA + data_file,python_dict = True)
    cluster_struct = list(data_dict.values())

    bolocam,err = clus_convert_bolocam(cluster_struct,verbose=verbose)
    if err:
        errmsg = str('clus_convert_bolocam exited with error: ' + err)
        if verbose:
            print(errmsg)
        return None, errmsg

    # Set the size of the image for later use
    naxis = bolocam[0]['deconvolved_image'][0].shape
    naxis = np.array(naxis)

    sziso = fits.PrimaryHDU()
    temphead = sziso.header

    # Set reference pizel position
    crpix1 = int(naxis[0] / 2)
    crpix2 = int(naxis[1] / 2)

    # Setting the tempheader for the bolocam data
    temphead.set('CRVAL1' , bolocam[0]['deconvolved_image_ra_j2000_deg'][0][crpix1,crpix2])
    temphead.set('CRVAL2' , bolocam[0]['deconvolved_image_dec_j2000_deg'][0][crpix1,crpix2])
    temphead.set('CRPIX1' , crpix1)
    temphead.set('CRPIX2' , crpix2)
    temphead.set('CD1_1' , -bolocam[0]['deconvolved_image_resolution_arcmin'][0] / 60.0)
    temphead.set('CD1_2' , 0)
    temphead.set('CD2_1' ,  0)
    temphead.set('CD2_2' , bolocam[0]['deconvolved_image_resolution_arcmin'][0] / 60.0)
    temphead.set('EPOCH' , 2000)
    temphead.set('EQUINOX' , 2000)
    temphead.set('CTYPE1' , 'RA---TAN')
    temphead.set('CTYPE2' , 'DEC--TAN')

    x = np.arange(naxis[0]) - 14.5
    y = np.arange(naxis[1]) - 14.5

    rad = np.zeros((naxis[0],naxis[1]))

    for ix in range(naxis[1]):
        rad[:,ix] = (np.tile(x[ix],naxis[0])**2+y**2)**(1/2)

    # Outer marks the indices of the outside of the circle used to add in the sz effect
    n = 0
    outer = []
    for i in range(len(rad)):
        for j in range(len(rad)):
            if rad[i,j] > 10:
                outer.append(n)
            n +=1

    # Set up the spectral shape of the sz effect to be appled to the 500um map
    if ib_test == False :
        # plt.imshow(bolocam[0]['deconvolved_image'][0])
        # plt.show()
        szmap = -1 * bolocam[0]['deconvolved_image'][0]
        szmap = (np.array(szmap)).flatten()
        szmap = szmap - np.mean(szmap[outer])
        szmap = [x/max(szmap) for x in szmap]
        # szmap = np.reshape(szmap, (naxis[0], naxis[1]))

    for imap in range(mapsize):
        # Applying the effect to the 500um band
        if imap == 2:
            yin_coeff = [2.50,1.91,2.26,3.99,1.36,2.42,1.59,1.90,3.99]
            yin = [x*1e-4 for x in yin_coeff]
            tin = [7.2,10.1,7.7,9.8,4.5,8.6,7.8,5.5,10.9]
            nu = 3e5 / clus_get_lambdas((maps[imap]['band']))
            dI,errmsg = clus_get_relsz(nu,y=yin,te=tin,vpec=0.0) # dI = [MJy/sr]
            print(dI / maps[imap]['calfac'])
            if errmsg:
                if verbose:
                    new_errmsg = 'Clus_get_relSZ exited with error'+errmsg
                return None, new_errmsg


            ''' This was necessary using old units for run_SZpack output '''
            # dI_converted = dI*((1.13*(maps[imap]['widtha'])/3600.)**2 * (pi/180)**2)

            if ib_test == True :
                szmap,err = IB_model(maps[imap],params,verbose)
                plt.imshow(szmap)
                plt.show()
                szmap = (np.array(szmap))
                naxis1 = szmap.shape[0]
                naxis2 = szmap.shape[1]
                szmap = szmap.flatten()
                szin = [x * dI / maps[imap]['calfac'] for x in szmap]
                szinp = np.reshape(szin,(naxis1,naxis2))

            else :
                # Combine the spectral shape of the SZ effect, and combine with the peak intensity
                # converted to Jy/pixel  ***Confirm these Units***
                szin = [x * dI  for x in szmap] #/ maps[imap]['calfac']
                szin = np.reshape(szin,(naxis[0],naxis[1]))

            # Have to interpolate to SPIRE map size
            # Requires fits objects to work
                hdu = fits.PrimaryHDU(szin,temphead)
                hdx = fits.PrimaryHDU(maps[imap]['signal'],maps[imap]['shead'])
                szinp = hcongrid(hdu.data,hdu.header,hdx.header)

                #Used to check the alligned sz effect image
                # sz = fits.PrimaryHDU(szinp,hdx.header)
                # sz.writeto('test.fits')

            # Combine the original signal with the sz effect
            '''TESTING : For now just passing through the sz signal only'''
            # maps[imap]['signal'] = maps[imap]['signal'] + szinp
            maps[imap]['signal'] = szinp
            # plt.imshow(szinp)
            # plt.show()

            # check the final output map
            # sz = fits.PrimaryHDU(maps[imap]['signal'],hdx.header)
            # sz.writeto('a0370_sim200PLW_sz_only.fits')

    return maps, None

if __name__ == '__main__' :
    yin_coeff = [2.50,1.91,2.26,3.99,1.36,2.42,1.59,1.90,3.99]
    yin = [x*1e-4 for x in yin_coeff]
    tin = [7.2,10.1,7.7,9.8,4.5,8.6,7.8,5.5,10.9]
    maps,err = clus_get_data('a0370')
    clus_add_sziso(maps,yin=yin,tin=tin)
