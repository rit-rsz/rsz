
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
# REVISION HISTORY : 3/6/2020 - VLB added in new Bolocam template for RXJ1347
################################################################################
import scipy.io
import numpy as np
import math
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
from IB_model import *
import matplotlib.pyplot as plt
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve_fft as convolve

def clus_add_sziso(maps,isim,yin,tin,params,
              verbose = 0, testflag=0, saveplot=0,nsim=0,clusname=None):

    errmsg = False

    # Now the bolocam data is in the SPIRE format
    # we can do a loop over the map adding in the false sz signal

    mapsize = len(maps)

    if verbose:
        print('Fetching BOLOCAM maps')

    ''' This is just until Jack gives us more new Bolocam templates... '''
    if clusname == 'rxj1347' :
        data_file = CLUSDATA + 'bolocam/new/' + maps[0]['name'] + '.fits'
        # cluster_struct = fits.getdata(data_file)
        hdul = fits.open(data_file)
        header = hdul[0].header
        naxis = hdul[0].data.shape
        # BCAMNORM in units of MJy/sr
        bolocam = clus_convert_bolocam(hdul[0].data,norm = header['BCAMNORM'],verbose=verbose,clusname=clusname)
    else :
        #fetch bolocam data from .sav files in bolocam directory.
        data_file = str('bolocam/data/' + maps[0]['name'] + '.sav')
        data_dict = scipy.io.readsav(CLUSDATA + data_file,python_dict = True)
        cluster_struct = list(data_dict.values())

        bolocam,err = clus_convert_bolocam(cluster_struct,verbose=verbose,clusname=clusname)
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

    if clusname == 'rxj1347' :
        # Setting the tempheader for the bolocam data
        temphead.set('CRVAL1' , header['CRVAL1']) # DEC in deg of ref pixel
        temphead.set('CRVAL2' , header['CRVAL2']) # RA in deg of ref pixel
        temphead.set('CRPIX1' , header['CRPIX1']) # ref pixel x
        temphead.set('CRPIX2' , header['CRPIX2']) # ref pixel y
        temphead.set('CD1_1' , header['CD1_1']) # Deg/pixel
        temphead.set('CD1_2' , header['CD1_2']) # Deg/pixel
        temphead.set('CD2_1' , header['CD2_1']) # Deg/pixel
        temphead.set('CD2_2' , header['CD2_2']) # Deg/pixel
        temphead.set('EPOCH' , 2000)
        temphead.set('EQUINOX' , 2000)
        temphead.set('CTYPE1' , header['CTYPE1']) # coord type
        temphead.set('CTYPE2' , header['CTYPE2']) # coord type

    else :
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

    # 14.5  = full size of the bolocam image in pixels ?
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
    if testflag == 0:
        if clusname == 'rxj1347' :
            szmap1 = hdul[0].data
        else :
            szmap1 = -1 * bolocam[0]['deconvolved_image'][0] # -1 is needed because SZ amp flips below 217 GHz
        szmap1 = (np.array(szmap1)).flatten()
        sz_mean = np.mean(szmap1[outer])
        szmap2 = [x - sz_mean for x in szmap1]
        sz_max = max(szmap2)
        szmap = [x/sz_max for x in szmap2]

    final_dI = []
    for imap in range(mapsize):
        # Applying the effect to the 500 um and 350 um bands.
        if imap == 2 or imap == 1:
            if testflag == 1:
                szmap,err = IB_model(maps[imap],params,verbose)
                szmap = np.array(szmap)
                plt.imshow(szmap)
                plt.colorbar().set_label(['Jy'])
                plt.title('Clus Add Sziso : IB Model for %s' %(maps[imap]['name']))
                plt.savefig(config.OUTPUT + 'add_sz/%s_ibmodel_%s_%s.png' %(maps[imap]['name'],maps[imap]['band'],isim))
                plt.clf()
                naxis = szmap.shape
                szmap = szmap.flatten()

            nu = 3e5 / clus_get_lambdas((maps[imap]['band']))

            if int(isim) == 0 : # dI only needs to be calculated once...
                dI,errmsg = clus_get_relsz(isim,nu,imap,y=yin,te=tin,vpec=0.0) # dI = [MJy/sr]
                np.save(config.OUTPUT + 'add_sz/sim_dI_%s.npy' %(maps[imap]['band']),dI)
            else :
                dI = np.load(config.OUTPUT + 'add_sz/sim_dI_%s.npy' %(maps[imap]['band']))

            if errmsg:
                if verbose:
                    new_errmsg = 'Clus_get_relSZ exited with error'+errmsg
                return None, new_errmsg

            # Combine the spectral shape of the SZ effect, and combine with the peak intensity
            # converted to Jy/beam  ***Confirm these Units***
            if int(isim) == 0 :
                szin = [x * dI / maps[imap]['calfac'] for x in szmap]
                final_dI.append(dI / maps[imap]['calfac'])
            else :
                szin = [x * dI / maps[imap]['calfac'] for x in szmap]
                final_dI.append(dI / maps[imap]['calfac'])

            szin = np.reshape(szin,(naxis[0],naxis[1]))

            # Have to interpolate to SPIRE map size
            hdu = fits.PrimaryHDU(szin,temphead)
            hdx = fits.PrimaryHDU(maps[imap]['signal'],maps[imap]['shead'])

            if testflag == 0 :
                szinp = hcongrid(hdu.data,hdu.header,hdx.header)
                # need to smooth output with SPIRE PSF
                fwhm = maps[imap]['widtha']
                pixscale = maps[imap]['pixsize']
                retext = round(fwhm * 5.0 / pixscale)
                if retext % 2 == 0:
                    retext += 1
                bmsigma = fwhm / math.sqrt(8 * math.log(2))
                beam = Gaussian2DKernel(bmsigma / pixscale, x_size=retext,y_size=retext, mode='oversample',factor=1)
                beam *= 1.0 / beam.array.max()
                out_map = convolve(szinp, beam, boundary='wrap')

                # Combine the original signal with the sz effect
                maps[imap]['signal'] = maps[imap]['signal'] + out_map
            else :
                # maps[imap]['signal'] = maps[imap]['signal'] + szin
                '''THIS IS JUST FOR TESTING WITH SZ SIGNAL'''
                maps[imap]['signal'] = maps[imap]['error'] + szin

            # Used to check the alligned sz effect image
            if saveplot:
                filename = config.OUTPUT + 'sim_sz/' + maps[imap]['name'] + '_sze_' + maps[imap]['band'] + '_' + str(isim) + '.png'
                plt.imshow(maps[imap]['signal'])
                plt.title('Clus Add Sziso : SZE + Signal Map %s' %(maps[imap]['band']))
                plt.colorbar().set_label('[Jy]')
                plt.savefig(filename)
                plt.clf()

    return maps, None, final_dI

if __name__ == '__main__' :
    yin_coeff = [2.50*1e-4]#,1.91,2.26,3.99,1.36,2.42,1.59,1.90,3.99]
    # yin = [x*1e-4 for x in yin_coeff]
    tin = [7.2]#,10.1,7.7,9.8,4.5,8.6,7.8,5.5,10.9]
    maps,err = clus_get_data('a0370')
    clus_add_sziso(maps,yin=yin,tin=tin)
