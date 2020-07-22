
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

def clus_add_sziso_new(maps,isim,yin=0,tin=0,params=None,
              verbose = 0, testflag=0, saveplot=0,clusname=None):

    #------------------------------------- This part is just for testing
    yin = 9.65*1e-4
    tin = 10.88

    errmsg = False

    # Now the bolocam data is in the SPIRE format
    # we can do a loop over the map adding in the false sz signal

    mapsize = len(maps)

    if verbose:
        print('Fetching BOLOCAM maps')

    ''' This is just until Jack gives us more new Bolocam templates... '''
    if clusname == 'rxj1347' :
        data_file = CLUSDATA + 'bolocam/new/' + maps[0]['name'] + '.fits'
        cluster_struct = fits.getdata(data_file)
        data = fits.open(data_file)
        header = data[0].header
        naxis = cluster_struct.shape
        # BCAMNORM in units of MJy/sr
        print(header['BCAMNORM'], 'bolo cam norm')
        bolocam,err = clus_convert_bolocam(cluster_struct,norm = header['BCAMNORM'],verbose=verbose,clusname=clusname)

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
        temphead = header

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
            szmap1 = bolocam
        else :
            szmap1 = -1 * bolocam[0]['deconvolved_image'][0] # -1 is needed because SZ amp flips below 217 GHz
            szmap1 = (np.array(szmap1)).flatten()
            sz_mean = np.mean(szmap1[outer])
            szmap2 = [x - sz_mean for x in szmap1]
            sz_max = max(szmap2)
            szmap = [x/sz_max for x in szmap2]

    final_dI = []
    if clusname == 'rxj1347' : # need to calculate the dI for Bolocam using new maps
        # if os.path.isfile(config.HOME + 'lookup/rxj1347_bol_lookup.npy') : # use lookup file instead of running szpack every time
        #     dI_bolo = np.load(config.HOME + 'lookup/rxj1347_bol_lookup.npy')
        # else :
        szmap1 = bolocam
        dI_bolo,thisx,JofXout,xout,errmsg = clus_get_relsz(isim,3e5 / clus_get_lambdas('BOLOCAM'),'BOLOCAM',y=yin,te=tin,vpec=0.0) # dI = [MJy/sr]
        # np.save(config.HOME + 'lookup/rxj1347_bol_lookup.npy',dI_bolo)
        szmap1 = np.array(szmap1).flatten()
        szmap1 = [x/dI_bolo for x in szmap1]

    for imap in range(mapsize):
        # Applying the effect to the 500 um and 350 um bands.
        # if imap == 2 or imap == 1:
        if testflag == 1:
            szmap,err = IB_model(maps[imap],params,verbose)
            szmap = np.array(szmap)
            plt.imshow(szmap,origin=0)
            plt.colorbar().set_label(['Jy'])
            plt.title('Clus Add Sziso : IB Model for %s' %(maps[imap]['name']))
            plt.savefig(config.OUTPUT + 'add_sz/%s_ibmodel_%s_%s.png' %(maps[imap]['name'],maps[imap]['band'],isim))
            plt.clf()
            naxis = szmap.shape
            szmap = szmap.flatten()

        nu = 3e5 / clus_get_lambdas((maps[imap]['band']))

        if int(isim) == 0 : # dI only needs to be calculated once...
            dI,thisx,JofXout,xout,errmsg = clus_get_relsz(isim,nu,imap,y=yin,te=tin,vpec=0.0) # dI = [MJy/sr]
            np.save(config.OUTPUT + 'add_sz/sim_dI_%s.npy' %(maps[imap]['band']), dI) #just the expected amplitude at spire wavelength
            np.save(config.OUTPUT + 'add_sz/input_dI_%s.npy' % (maps[imap]['band']), dI * -.94 / dI_bolo) #input amplitude saved for comparison in fit_sz
        else :
            dI = np.load(config.OUTPUT + 'add_sz/sim_dI_%s.npy' %(maps[imap]['band']))

        if errmsg:
            if verbose:
                new_errmsg = 'Clus_get_relSZ exited with error'+errmsg
            return None, new_errmsg

        # Combine the spectral shape of the SZ effect, and combine with the peak intensity
        # converted to Jy/beam

        hdtemp = fits.PrimaryHDU(szmap1, temphead)
        hdtemp.writeto(config.OUTPUT + 'sz_fit/bolo_fit_template%s.fits' % isim, overwrite=True)

        szin = [x * dI/ maps[imap]['calfac'] for x in szmap1]  #
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

            '''THIS IS JUST FOR TESTING WITH SZ SIGNAL'''
            # maps[imap]['signal'] = out_map
        elif testflag == 2 :

            szinp = hcongrid(hdu.data,hdu.header,hdx.header)
            # need to smooth output with SPIRE PSF
            fwhm = maps[imap]['widtha']
            pixscale = maps[imap]['pixsize']
            retext = round(fwhm * 5.0 / pixscale)
            if retext % 2 == 0:
                retext += 1
            bmsigma = fwhm / math.sqrt(8 * math.log(2))
            beam = np.asarray(Gaussian2DKernel(bmsigma / pixscale, x_size=retext,y_size=retext, mode='oversample',factor=1))
            n_beam = beam / np.sum(beam) #the input bolocam map is a surface brightness so we need to normalize over area
            out_map = convolve(szinp, n_beam, boundary='wrap')

            maps[imap]['signal'] = maps[imap]['signal'] + out_map
            # maps[imap]['srcrm'] = out_map #for testing
            # maps[imap]['signal'] = maps[imap]['error'] + out_map

        # Used to check the alligned sz effect image
        if saveplot:
            filename = config.OUTPUT + 'sim_sz/' + maps[imap]['name'] + '_sze+signal_' + maps[imap]['band'] + '_' + str(isim) + '.png'
            plt.imshow(maps[imap]['signal'],origin=0)
            plt.title('Clus Add Sziso : SZE + Signal Map %s' %(maps[imap]['band']))
            plt.colorbar().set_label('[Jy]')
            plt.clim(-0.1,0.1)
            plt.savefig(filename)
            plt.clf()

            filename = config.OUTPUT + 'sim_sz/' + maps[imap]['name'] + '_sze_' + maps[imap]['band'] + '_' + str(isim) + '.png'
            plt.imshow(out_map, origin=0)
            plt.title('Clus Add Sziso : SZE + Signal Map %s' %(maps[imap]['band']))
            plt.clim(-0.1,0.1)
            plt.colorbar().set_label('[Jy]')
            plt.savefig(filename)
            plt.clf()

            savefile = config.OUTPUT + 'sim_sz/' + maps[imap]['name'] + '_sze_' + maps[imap]['band'] + '_' + str(isim) + '.fits'
            hda = fits.PrimaryHDU(maps[imap]['signal'],maps[imap]['shead'])
            hda.writeto(savefile,overwrite=True)

    return maps, None, final_dI

if __name__ == '__main__' :
    sys.path.append('../source_handling')
    from clus_get_data import *
    from clus_get_clusparams import *
    # yin = 14.08*1e-4
    yin = 9.65*1e-4
    tin = 10.88
    maps,err = clus_get_data('rxj1347','0',sgen=3)
    params, err = clus_get_clusparams('rxj1347',0,verbose=0)
    clus_add_sziso_new(maps,0,yin=yin,tin=tin,params=params,verbose = 0, testflag=0, saveplot=1,clusname='rxj1347')
