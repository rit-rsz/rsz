import scipy.io
import numpy as np
import math, os, sys
sys.path.append('../utilities')
from config import *
sys.path.append('../sz')
sys.path.append('../source_handling')
from astropy.io import fits
from FITS_tools.hcongrid import hcongrid
import matplotlib.pyplot as plt


def bolocam_mask(maps) :

    #fetch bolocam data from .sav files in bolocam directory.
    data_file = str('bolocam/data/' + maps[0]['name'] + '.sav')
    data_dict = scipy.io.readsav(CLUSDATA + data_file,python_dict = True)
    bolocam = list(data_dict.values())
    # data = bolocam[0]['deconvolved_image']
    naxis = bolocam[0]['deconvolved_image'][0].shape
    naxis = np.array(naxis)
    sziso = fits.PrimaryHDU()
    temphead = sziso.header
    crpix1 = int(naxis[0] / 2.0)
    crpix2 = int(naxis[1] / 2.0)

    bolo_mask = np.ones((naxis[0],naxis[1]))

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

    for col in range(3):
        # retrieve error map from sim
        mask = maps[col]['error']
        err_mask = np.zeros((mask.shape[0],mask.shape[1]))

        # scale bolocam mask up to SPIRE pixel sizes
        hdu = fits.PrimaryHDU(bolo_mask,temphead)
        hdx = fits.PrimaryHDU(maps[col]['signal'],maps[col]['shead'])
        temp_mask = hcongrid(hdu.data,hdu.header,hdx.header)
        final_mask = err_mask + temp_mask
        # plt.imshow(final_mask,origin=0)
        # plt.colorbar()
        # plt.savefig('temp_bolo_mask_%s.png'%(col))
        # plt.clf()

        # make all pixels outside of bolo_mask zeros
        for j in range(mask.shape[0]):
            for k in range(mask.shape[1]):
                if final_mask[j,k] == 0 :
                    mask[j,k] = 0
                else :
                    mask[j,k] = 1

        hda = fits.PrimaryHDU(mask,maps[col]['shead'])
        hda.writeto('bolocam_mask_%s.fits'%(maps[col]['band']),overwrite=True)
        plt.imshow(mask,origin=0)
        plt.colorbar()
        plt.savefig('final_bolo_mask_%s.png'%(col))
        plt.clf()
