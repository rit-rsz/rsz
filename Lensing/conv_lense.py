# python
import numpy as np
import matplotlib.pyplot as plt
import math
import sys, time
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve_fft as convolve

def conv_lense():
    psw = 'lense_template_PSW.fits'
    pmw = 'lense_template_PMW.fits'
    plw = 'lense_template_PLW.fits'
    psw_hdul = fits.open(psw)
    psw_header = psw_hdul[0].header
    psw_data = psw_hdul[0].data
    pmw_hdul = fits.open(pmw)
    pmw_header = pmw_hdul[0].header
    pmw_data = pmw_hdul[0].data
    plw_hdul = fits.open(plw)
    plw_header = plw_hdul[0].header
    plw_data = plw_hdul[0].data
    fwhm = [17.6, 23.9, 35.2]
    pixsize = [6.0, 8.33333, 12.0]
    bands = ['PSW','PMW','PLW']
    conv_map(psw_data,fwhm[0],pixsize[0],bands[0],psw_header)
    conv_map(pmw_data,fwhm[1],pixsize[1],bands[1],pmw_header)
    conv_map(plw_data,fwhm[2],pixsize[2],bands[2],plw_header)
    return

def conv_map(in_map,fwhm,pixsize,bands,header):
    beam = get_gauss_beam(fwhm,pixsize,oversamp=5)
    outmap = convolve(in_map, beam, boundary='wrap')
    plt.imshow(outmap,origin=0)
    plt.clim([-0.02,0.02])
    plt.colorbar()
    plt.title('Beam Convolved Lense Temp %s'%(bands))
    plt.savefig('conv_lense_template_%s.png'%(bands))
    plt.clf()
    hda = fits.PrimaryHDU(outmap,header)
    hda.writeto('conv_lense_template_%s.fits'%(bands),overwrite=True)
    return


def get_gauss_beam(fwhm, pixscale, nfwhm=5.0, oversamp=1):
    retext = round(fwhm * nfwhm / pixscale)
    if retext % 2 == 0:
        retext += 1

    bmsigma = fwhm / math.sqrt(8 * math.log(2))

    beam = Gaussian2DKernel(bmsigma / pixscale, x_size=retext,
                            y_size=retext, mode='oversample',
                            factor=oversamp)
    beam *= 1.0 / beam.array.max()
    return beam

if __name__ == '__main__' :
    conv_lense()
