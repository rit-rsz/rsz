import numpy as np
import matplotlib.pyplot as plt
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

maps,err = clus_get_data('rxj1347',0)
bands = ['PSW','PMW','PLW']

# for col in range(3):
#     mapsize = maps[col]['error'].shape
#     img = maps[col]['error']
#     stddev = 24 / 2.355
#     kernel = makeGaussian(8, 8, fwhm= 4)
#     kernel_full_size = padGaussian(img, kernel)
#
#     #converting these to lists because the error message is annoying.
#     img = img.tolist()
#     kernel_full_size = kernel_full_size.tolist()
#     fixed_image = convolve_fft(img, kernel_full_size)
#     plt.imshow(fixed_image,origin=0)
#     plt.clim([-0.005,0.005])
#     plt.colorbar()
#     plt.title('Gaussian Smoothed SPIRE Error %s'%(bands[col]))
#     plt.savefig('convolve_error_%s.png'%(bands[col]))
#     plt.clf()
# exit()


mask_psw = np.array(clus_make_noise_mask(maps,0))
mask_pmw = clus_make_noise_mask(maps,1)
# mask_pmw = np.array([x*10 for x in mask_pmw])
mask_plw = clus_make_noise_mask(maps,2)
# mask_plw = np.array([x*100 for x in mask_plw])

for j in range(maps[0]['signal'].shape[0]):
    for k in range(maps[0]['signal'].shape[1]):
        if mask_psw[j,k] == 1:
            maps[0]['error'][j,k] = np.nan
        else :
            maps[0]['error'][j,k] = mask_psw[j,k]

for j in range(maps[1]['signal'].shape[0]):
    for k in range(maps[1]['signal'].shape[1]):
        if mask_pmw[j,k] == 1:
            maps[1]['error'][j,k] = np.nan
        else :
            maps[1]['error'][j,k] = mask_pmw[j,k]

for j in range(maps[2]['signal'].shape[0]):
    for k in range(maps[2]['signal'].shape[1]):
        if mask_plw[j,k] == 1:
            maps[2]['error'][j,k] = np.nan
        else :
            maps[2]['error'][j,k] = mask_plw[j,k]

bolo_pmw = hcongrid(maps[1]['error'],maps[1]['shead'],maps[2]['shead'])
bolo_plw = hcongrid(maps[0]['error'],maps[0]['shead'],maps[2]['shead'])

comp_plw = bolo_pmw + bolo_plw

plt.imshow(comp_plw,origin=0)
# plt.imshow(comp,origin=0,norm=LogNorm(vmin=1, vmax=100))
plt.colorbar()
# plt.clim([-0.005,0.005])
plt.title('Compilation Mask')
plt.savefig('comp_mask_final.png')
plt.clf()

# set final PLW maskt to whatever is non zero here
for j in range(maps[2]['signal'].shape[0]):
    for k in range(maps[2]['signal'].shape[1]):
        if maps[0]['error'][j,k] == np.nan or maps[1]['error'][j,k] == np.nan or maps[2]['error'][j,k] == np.nan:
            maps[2]['error'][j,k] = np.nan
        else :
            maps[2]['error'][j,k] = 1

plt.imshow(comp_plw,origin=0)
# plt.imshow(comp,origin=0,norm=LogNorm(vmin=1, vmax=100))
plt.colorbar()
# plt.clim([-0.005,0.005])
plt.title('Compilation Mask Inverted')
plt.savefig('comp_mask.png')
plt.clf()

''' ###############################################################
# bolocam = fits.open(config.HOME + 'filtered_image.fits')
# header = bolocam[0].header
# temphead = fits.PrimaryHDU()
# bolohead= temphead.header
# # Setting the tempheader for the bolocam data
# bolohead.set('CRVAL1' , header['CRVAL1']) # DEC in deg of ref pixel
# bolohead.set('CRVAL2' , header['CRVAL2']) # RA in deg of ref pixel
# bolohead.set('CRPIX1' , header['CRPIX1']) # ref pixel x
# bolohead.set('CRPIX2' , header['CRPIX2']) # ref pixel y
# bolohead.set('CD1_1' , header['CD1_1']) # Deg/pixel
# bolohead.set('CD1_2' , header['CD1_2']) # Deg/pixel
# bolohead.set('CD2_1' , header['CD2_1']) # Deg/pixel
# bolohead.set('CD2_2' , header['CD2_2']) # Deg/pixel
# bolohead.set('EPOCH' , 2000)
# bolohead.set('EQUINOX' , 2000)
# bolohead.set('CTYPE1' , header['CTYPE1']) # coord type
# bolohead.set('CTYPE2' , header['CTYPE2']) # coord type

# bolo_shape = bolocam[0].data.shape
# bolo = np.ones((bolo_shape[0],bolo_shape[1]))
# hdu = fits.PrimaryHDU(bolo,bolohead)
#
# hdx_psw = fits.PrimaryHDU(mask_psw,maps[0]['shead'])
# bolo_psw = hcongrid(hdu.data,hdu.header,hdx_psw.header)
#
# hdx_pmw = fits.PrimaryHDU(mask_pmw,maps[1]['shead'])
# bolo_pmw = hcongrid(hdu.data,hdu.header,hdx_pmw.header)
#
# hdx_plw = fits.PrimaryHDU(mask_plw,maps[2]['shead'])
# bolo_plw = hcongrid(hdu.data,hdu.header,hdx_plw.header)

# final_psw = bolo_psw * mask_psw
# final_pmw = bolo_pmw * mask_pmw
# final_plw = bolo_plw * mask_plw

# hda = fits.PrimaryHDU(final_psw,maps[0]['shead'])
# hda.writeto('bolocam_mask_psw.fits',overwrite=True)
# hdb = fits.PrimaryHDU(final_pmw,maps[1]['shead'])
# hdb.writeto('bolocam_mask_pmw.fits',overwrite=True)
# hdc = fits.PrimaryHDU(final_plw,maps[2]['shead'])
# hdc.writeto('bolocam_mask_plw.fits',overwrite=True)

# plt.imshow(final_psw,origin=0)
# plt.colorbar()
# # plt.clim([-0.005,0.005])
# plt.title('BOLOCAM + NOISE MASK PSW')
# plt.savefig('bolocam_mask_psw.png')
# plt.clf()
#
#
# plt.imshow(final_pmw,origin=0)
# plt.colorbar()
# # plt.clim([-0.005,0.005])
# plt.title('BOLOCAM + NOISE MASK PMW')
# plt.savefig('bolocam_mask_pmw.png')
# plt.clf()
#
# plt.imshow(final_plw,origin=0)
# plt.colorbar()
# # plt.clim([-0.005,0.005])
# plt.title('BOLOCAM + NOISE MASK PLW')
# plt.savefig('bolocam_mask_plw.png')
# plt.clf()
###############################################################   '''
