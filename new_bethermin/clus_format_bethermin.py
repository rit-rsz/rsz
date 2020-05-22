############################################################################
#
# NAME : clus_format_bethermin.py
# DATE : October 3, 2019
# AUTHOR : Victoria Butler
# PURPOSE : takes the catalog output by Alex Conley's implementation of
#           the bethermin model and puts it into the format that lenstool wants.
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
# OUTPUTS :
# REVISION HISTORY :

############################################################################
import sys, os, math
sys.path.append('../utilities')
import config
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from FITS_tools.hcongrid import hcongrid
from gaussian import makeGaussian
import math
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve_fft as convolve
from astropy.wcs.utils import skycoord_to_pixel
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
sys.path.append('../multiband_pcat')
from image_eval import psf_poly_fit, image_model_eval

def clus_format_bethermin(icol,sim_map,maps,map_size,band,clusname,pixsize,fwhm,\
                          fluxcut=0,zzero=0,superplot=0,savemaps=0,genbethermin=0):

    # trimming num sources for lenstool limit
    msrc = 50000 - 1

    if genbethermin :
        refx = maps[0]['shead']['CRPIX1']
        refy = maps[0]['shead']['CRPIX2'] # needs to be based on PSW map since all coordinates are referenced by pixel to maps[0] size
        # 3,4,5 are 250,350,500 truthtables
        cat = sim_map[icol+3]
        nsrc = len(cat['fluxdens']) # pulls len of truthtable

        # massage data into new arrays
        xpos = sim_map[icol+3]['x']
        ypos = sim_map[icol+3]['y']
        zpos = sim_map[icol+3]['z']
        outflux = sim_map[-1]['fluxdens'][:,icol] # all three tables are the same, so just use the last one
        outx = [pixsize * (x - refx) for x in xpos]
        outy = [pixsize * (y - refy) for y in ypos]

    else :
        refx = maps[icol]['shead']['CRPIX1']
        refy = maps[icol]['shead']['CRPIX2']
        ra = sim_map[icol].item().get('RA')
        dec = sim_map[icol].item().get('DEC')
        outflux = sim_map[icol].item().get('Flux')
        zpos = sim_map[icol].item().get('Redshift')
        nsrc = len(zpos)
        x_pos = maps[icol]['signal'].shape[1]
        y_pos = maps[icol]['signal'].shape[0]
        min_ra = min(ra)
        min_dec = min(dec)
        outx = np.array([(((i-min_ra) * 3600.0 ) - (refx*pixsize)) for i in ra])
        outy = np.array([(((j-min_dec) * 3600.0 ) - (refy*pixsize)) for j in dec])

    if superplot :
        plt.scatter(outx,outy,s=2,c=outflux)
        plt.colorbar()
        plt.title('Bethermin SIM (pre-format)')
        plt.show()
        plt.clf()

    outz = [float(np.ceil(10.0 * z)) / 10.0 for z in zpos]

    # lets do some fluxcuts
    if fluxcut > 0.0 :
        print('Cutting input catalog at flux %s' %(fluxcut))
        for i in range(len(outflux)):
            if outflux[i] <= fluxcut :
                np.delete(outflux,i)
                np.delete(outx,i)
                np.delete(outy,i)
                np.delete(outz,i)
        nsrc = len(outflux)

    savex = []
    savey = []
    savez = []
    savef = []
    index = []
    if zzero > 0 :
        for j in range(len(outflux)):
            if outz[j] <= zzero :
                index.append(j)
                savef.append(outflux[j])
                savex.append(outx[j])
                savey.append(outy[j])
                savez.append(outz[j])
        retcat = {'x':savex,'y':savey,'z':savez,'f':savef}
        nsrc = len(outflux)
        coutx = np.delete(np.asarray(outx),index)
        couty = np.delete(np.asarray(outy),index)
        coutz = np.delete(np.asarray(outz),index)
        coutflux = np.delete(np.asarray(outflux),index)

    # sort according to brightness due to lenstool limitations
    # lambda pair: pair[0] tells sorted to use outflux as the sorting key
    outx = [x for _,x in sorted(zip(coutflux,coutx), key = lambda pair: pair[0], reverse=True)]
    outy = [y for _,y in sorted(zip(coutflux,couty), key = lambda pair: pair[0], reverse=True)]
    outz = [z for _,z in sorted(zip(coutflux,coutz), key = lambda pair: pair[0], reverse=True)]
    outflux = sorted(coutflux, reverse=True)

    if superplot :
        plt.scatter(outx,outy,s=2,c=outflux)
        plt.colorbar()
        plt.title('Bethermin SIM (pre-format)')
        plt.show()
        plt.clf()

    # truncate to the msrc brightest sources
    if msrc < nsrc :
        outflux = outflux[0:msrc]
        outx = outx[0:msrc]
        outy = outy[0:msrc]
        outz = outz[0:msrc]

    # now sort according to z
    coutflux = [f for _,f in sorted(zip(outz,outflux), key = lambda pair: pair[0])]
    coutx = [x for _,x in sorted(zip(outz,outx), key = lambda pair: pair[0])]
    couty = [y for _,y in sorted(zip(outz,outy), key = lambda pair: pair[0])]
    coutz = sorted(outz)

    if savemaps:
        orig_length = len(coutx)
        mapy = maps[icol]['signal'].shape[0]
        mapx = maps[icol]['signal'].shape[1]
        x = [(i/pixsize) + refx for i in coutx]
        y = [(j/pixsize) + refy for j in couty]

        x = np.array(x,dtype=np.float32)
        y = np.array(y,dtype=np.float32)
        flux = np.array(coutflux,dtype=np.float32)
        psf, cf, nc, nbin = get_gaussian_psf_template(fwhm,pixel_fwhm=3.) # assumes pixel fwhm is 3 pixels in each band
        if x_pos > y_pos :
            sim_map = image_model_eval(x, y, nc*flux, 0.0, (x_pos,x_pos), int(nc), cf)
        else :
            sim_map = image_model_eval(x, y, nc*flux, 0.0, (y_pos,y_pos), int(nc), cf)
        plt.imshow(sim_map,origin=0)
        plt.savefig(config.SIM + 'nonlensedmap_' + clusname + '_' + band + '.png')
        plt.colorbar()
        plt.clf()
        hdx = fits.PrimaryHDU(maps[icol]['signal'],maps[icol]['shead'])
        sz = fits.PrimaryHDU(sim_map,hdx.header)
        sz.writeto(config.SIMBOX + 'nonlensedmap_' + clusname + '_' + band + '.fits',overwrite=True)

    # magnitude instead of flux in Jy
    outmag = [-2.5 * np.log10(x) for x in coutflux]
    # if superplot and sim_map != None:
        plt.imshow(sim_map,extent=(0,300,0,300),clim=[0.0,0.15],origin=0)
        plt.colorbar()
        plt.title('clus_format_bethermin: Non-lensed map')
        plt.show()
        plt.clf()
    # else :
    #     plt.scatter(houtx,houty,s=2,c=houtflux)
    #     plt.colorbar()
    #     plt.title('end of format bethermin')
    #     plt.show()

    #write everything to file for lenstool to ingest
    lensfile = (config.HOME + 'model/' + clusname + '/' + clusname + '_cat.cat')
    with open(lensfile,'w') as f :
        f.write('#REFERENCE 3 %.6f %.6f \n' %(maps[icol]['shead']['CRVAL1'], maps[icol]['shead']['CRVAL2']))
        for k in range(len(outmag)):
            f.write('%i %.3f %.3f 0.5 0.5 0.0 %0.6f %0.6f \n' \
                    %(k,coutx[k],couty[k],coutz[k],outmag[k]))
        f.close()

    truthtable = {'x': coutx, 'y': couty,
                  'z': coutz, 'mag': outmag}

    return retcat, truthtable, sim_map

def get_gaussian_psf_template(fwhm,pixel_fwhm=3., nbin=5):
    nc = nbin**2
    psfnew = Gaussian2DKernel(pixel_fwhm/2.355*nbin, x_size=125, y_size=125).array.astype(np.float32)
    # psfnew2 = psfnew / np.max(psfnew)  * nc
    cf = psf_poly_fit(psfnew, nbin=nbin)
    return psfnew, cf, nc, nbin
