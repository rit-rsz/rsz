################################################################################
<<<<<<< Updated upstream
# NAME : clus_make_bolo_mask
# DATE STARTED : May 8, 2020
# AUTHORS : Victoria Butler, Benjamin Vaughan
# PURPOSE : This script interpolates IRAS coordinates to SPIRE coordinates
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
# OUTPUTS :
# REVISION HISTORY :
################################################################################import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from math import *
import sys, time
sys.path.append('source_handling')
sys.path.append('utilities')
from clus_get_data import clus_get_data
from astropy.io import fits
from gaussian import makeGaussian
from astropy.convolution import convolve_fft
from gaussian import padGaussian
from FITS_tools.hcongrid import hcongrid
from clus_make_noise_mask import clus_make_noise_mask
import config
from matplotlib.colors import LogNorm
from gaussian import makeGaussian
from astropy.convolution import convolve_fft
from gaussian import padGaussian
from clus_interp_map import interp_band_to_band


def mask_model(maps):
    #PIXSIZE -> LArge to SMALL PSW, PMW, PLW
    for i in range(len(maps)-1):
        base_mask   = noise_mask(maps, i)
        maps[i]['mask'] = np.round(interp_band_to_band(base_mask, maps[2], maps[i]))

    temp_plw_mask = noise_mask(maps, 2)


    composite = temp_plw_mask * maps[1]['mask'] * maps[0]['mask']


    for i in range(len(maps)):
        maps[i]['mask'] = np.round(interp_band_to_band(composite, maps[i], maps[2]))

    for i in range(len(maps)):
        plt.imshow(maps[i]['mask'])
        plt.savefig('round_test_mask%s' % maps[i]['band'])
        plt.clf()

def bolocam_mask():
    bolocam = fits.open(config.OUTPUT + 'bolocam_master_template.fits')
    header = bolocam[0].header
    temphead = fits.PrimaryHDU()
    bolohead = temphead.header
    # Setting the tempheader for the bolocam data
    bolohead.set('CRVAL1' , header['CRVAL1']) # DEC in deg of ref pixel
    bolohead.set('CRVAL2' , header['CRVAL2']) # RA in deg of ref pixel
    bolohead.set('CRPIX1' , header['CRPIX1']) # ref pixel x
    bolohead.set('CRPIX2' , header['CRPIX2']) # ref pixel y
    bolohead.set('CD1_1' , header['CD1_1']) # Deg/pixel
    bolohead.set('CD1_2' , header['CD1_2']) # Deg/pixel
    bolohead.set('CD2_1' , header['CD2_1']) # Deg/pixel
    bolohead.set('CD2_2' , header['CD2_2']) # Deg/pixel
    bolohead.set('EPOCH' , 2000)
    bolohead.set('EQUINOX' , 2000)
    bolohead.set('CTYPE1' , header['CTYPE1']) # coord type
    bolohead.set('CTYPE2' , header['CTYPE2']) # coord type

    bolo_shape = bolocam[0].data.shape
    bolo = np.ones((bolo_shape[0],bolo_shape[1]))
    hdu = fits.PrimaryHDU(bolo,bolohead)
    bolo_psw = hcongrid(hdu.data,hdu.header,maps[0]['shead'])
    bolo_pmw = hcongrid(hdu.data,hdu.header,maps[1]['shead'])
    bolo_plw = hcongrid(hdu.data,hdu.header,maps[2]['shead'])

    # plt.imshow(bolo_plw,origin=0)
    # plt.clim([0,1])
    # plt.colorbar()
    # plt.savefig('mystery_mask.png')
    # exit()


def noise_mask(maps,col):
    mapsize = maps[col]['error'].shape
    img = maps[col]['error']
    stddev = 24 / 2.355
    kernel = makeGaussian(8, 8, fwhm= 4)
    kernel_full_size = padGaussian(img, kernel)

    #converting these to lists because the error message is annoying.
    img = img.tolist()
    kernel_full_size = kernel_full_size.tolist()
    fixed_image = convolve_fft(img, kernel_full_size)
    noise_thresh = 5.0 * fixed_image[floor(mapsize[0]/2),floor(mapsize[1]/2)]
    print('noise threshold :',noise_thresh)

    mask = np.zeros(mapsize, dtype=int)
    for j in range(mapsize[0]):
        for k in range(mapsize[1]):
            if fixed_image[j,k] == 0 or fixed_image[j,k] <= 1e-10 or np.isfinite(fixed_image[j,k]) == False:
                fixed_image[j,k] = noise_thresh + 1e-3
            if fixed_image[j,k] <= noise_thresh:
                mask[j,k] = 1

    ''' Temp fix for edge effects '''
    mask[0:10,:] = 0 # fix top 10 rows
    mask[-10:mapsize[0],:] = 0  # fix bottom 10 rows
    mask[:,0:10] = 0 # fix left 10 rows
    mask[:,-10:mapsize[1]] = 0 # fix right 10 rows
    return np.asarray(mask)


if __name__ == '__main__':
    maps,err = clus_get_data('rxj1347',0)
    bands = ['PSW','PMW','PLW']
    mask_model(maps)


# final_psw = (bolo_psw * maps[0]['error']).flatten()
# final_pmw = (bolo_pmw * maps[1]['error']).flatten()
# final_plw = (bolo_plw * maps[2]['error']).flatten()
#
# # we have to do this yet again, because hcongrid leaves weird values that aren't 1 or 0
# new_psw = np.array([round(x) for x in final_psw]).reshape(maps[0]['signal'].shape[0],maps[0]['signal'].shape[1])
# new_pmw = np.array([round(y) for y in final_pmw]).reshape(maps[1]['signal'].shape[0],maps[1]['signal'].shape[1])
# new_plw = np.array([round(z) for z in final_plw]).reshape(maps[2]['signal'].shape[0],maps[2]['signal'].shape[1])
#
# ''' Save bolocam masked object to sim files '''


# for i in range(100):
#     sim_maps,err = clus_get_data('rxj1347', i , sgen = 3)
#     sim_maps[0]['mask'] = bolo_psw
#     sim_maps[1]['mask'] = bolo_pmw
#     sim_maps[2]['mask'] = bolo_plw
#
#     for j in range(3) :
#         savefile = config.CLUSSIMS + 'rxj1347/rxj1347_' + bands[j] + '_sim03%02d.fits' %(i)
#         hda = fits.PrimaryHDU(sim_maps[j]['signal'],sim_maps[j]['shead'])
#         hdul = fits.HDUList(hdus=hda)
#         hdn = fits.ImageHDU(sim_maps[j]['error'],sim_maps[j]['shead'])
#         hde = fits.ImageHDU(sim_maps[j]['exp'].data,sim_maps[j]['shead'])
#         hdm = fits.ImageHDU(sim_maps[j]['mask'],sim_maps[j]['shead'])
#         hdul.append(hdn)
#         hdul.append(hde)
#         hdul.append(hdm)
#         hdul.writeto(savefile,overwrite=True)
''' ################################################ '''

''' ################ Output images ################'''
hda = fits.PrimaryHDU(bolo_psw,maps[0]['shead'])
hda.writeto('bolocam_mask_PSW.fits',overwrite=True)
hdb = fits.PrimaryHDU(bolo_pmw,maps[1]['shead'])
hdb.writeto('bolocam_mask_PMW.fits',overwrite=True)
hdc = fits.PrimaryHDU(bolo_plw,maps[2]['shead'])
hdc.writeto('bolocam_mask_PLW.fits',overwrite=True)

# plt.contour(bolo_psw, 1, colors='red')
# plt.imshow(maps[0]['error'],origin=0)
# plt.colorbar()
# plt.clim([0,1])
# plt.title('BOLOCAM + NOISE MASK PSW')
# plt.savefig('bolocam_mask_psw.png')
# plt.clf()
#
# plt.contour(bolo_pmw, 1, colors='red')
# plt.imshow(maps[1]['error'],origin=0)
# plt.colorbar()
# plt.clim([0,1])
# plt.title('BOLOCAM + NOISE MASK PMW')
# plt.savefig('bolocam_mask_pmw.png')
# plt.clf()
#
# plt.contour(bolo_plw, 1, colors='red')
# plt.imshow(maps[2]['error'],origin=0)
# plt.colorbar()
# plt.clim([0,1])
# plt.title('BOLOCAM + NOISE MASK PLW')
# plt.savefig('bolocam_mask_plw.png')
# plt.clf()

# new_maps,err = clus_get_data('rxj1347',0)
#
# hda = fits.PrimaryHDU(new_maps[0]['signal'],new_maps[0]['shead'])
# hdul = fits.HDUList(hdus=hda)
# hdn = fits.ImageHDU(new_maps[0]['error'],new_maps[0]['shead'])
# hdm = fits.ImageHDU(new_psw,new_maps[0]['shead'])
# hdul.append(hdn)
# hdul.append(hdm)
# hdul.writeto('final_comp_psw.fits',overwrite=True)
#
# hda = fits.PrimaryHDU(new_maps[1]['signal'],new_maps[1]['shead'])
# hdul = fits.HDUList(hdus=hda)
# hdn = fits.ImageHDU(new_maps[1]['error'],new_maps[1]['shead'])
# hdm = fits.ImageHDU(new_pmw,new_maps[1]['shead'])
# hdul.append(hdn)
# hdul.append(hdm)
# hdul.writeto('final_comp_pmw.fits',overwrite=True)
#
# hda = fits.PrimaryHDU(new_maps[2]['signal'],new_maps[2]['shead'])
# hdul = fits.HDUList(hdus=hda)
# hdn = fits.ImageHDU(new_maps[2]['error'],new_maps[2]['shead'])
# hdm = fits.ImageHDU(new_plw,new_maps[2]['shead'])
# hdul.append(hdn)
# hdul.append(hdm)
# hdul.writeto('final_comp_plw.fits',overwrite=True)﻿
''' ####################################################################### '''
