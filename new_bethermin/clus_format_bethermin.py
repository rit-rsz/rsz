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
        x_pos = math.floor(maps[icol]['signal'].shape[0])
        y_pos = math.floor(maps[icol]['signal'].shape[1])
        min_ra = min(ra)
        min_dec = min(dec)
        # trim sources outside of the real SPIRE map size
        rem = []
        for k in range(len(ra)):
            if ((ra[k]-min_ra) * 3600.0 / pixsize) > x_pos or ((dec[k]-min_dec) * 3600.0 / pixsize) > y_pos :
                rem.append(k)
        ra = np.delete(ra,rem)
        dec = np.delete(dec,rem)
        outflux = np.delete(outflux,rem)
        zpos = np.delete(zpos,rem)
        # outx = np.array([((i-min_ra) * 3600.0 / pixsize) for i in ra])
        # outy = np.array([((j-min_dec) * 3600.0 / pixsize) for j in dec])
        outx = np.array([((i) * 3600.0 / pixsize)-(x_pos/2.0) for i in ra])
        outy = np.array([((j) * 3600.0 / pixsize)-(y_pos/2.0) for j in dec])

    # if superplot :
    plt.scatter(outx,outy,s=2,c=outflux)
    plt.colorbar()
    plt.title('Bethermin SIM (pre-format)')
    plt.savefig('format_bethermin_%s.png' %(band))
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
    outflux = sorted(outflux, reverse=True)

    # truncate to the msrc brightest sources
    if msrc < nsrc :
        outflux = outflux[0:msrc]
        outx = outx[0:msrc]
        outy = outy[0:msrc]
        outz = outz[0:msrc]

    # now sort according to z
    houtflux = [f for _,f in sorted(zip(outz,outflux), key = lambda pair: pair[0])]
    houtx = [x for _,x in sorted(zip(outz,outx), key = lambda pair: pair[0])]
    houty = [y for _,y in sorted(zip(outz,outy), key = lambda pair: pair[0])]
    houtz = sorted(outz)

    if savemaps:
        cmap = np.zeros((map_size,map_size),dtype=np.float32)
        xf = np.floor(houtx)
        yf = np.floor(houty)
        nx, ny = cmap.shape
        np.place(xf, xf > nx-1, nx-1)
        np.place(yf, yf > ny-1, ny-1)
        for cx, cy, cf in zip(xf, yf, houtflux):
            cmap[int(cy), int(cx)] += cf  # Note transpose

        beam = get_gauss_beam(fwhm,pixsize,band,oversamp=5)
        sim_map = convolve(cmap, beam, boundary='wrap')

        hdx = fits.PrimaryHDU(maps[icol]['signal'],maps[icol]['shead'])
        sz = fits.PrimaryHDU(sim_map,hdx.header)
        sz.writeto(config.SIMBOX + 'nonlensedmap_' + clusname + '_' + band + '.fits',overwrite=True)

    # magnitude instead of flux in Jy
    outmag = [-2.5 * np.log10(x) for x in houtflux]

    if superplot and sim_map != None:
        plt.imshow(sim_map,extent=(0,300,0,300),clim=[0.0,0.15],origin=0)
        plt.colorbar()
        plt.title('clus_format_bethermin: Non-lensed map')
        plt.show()
    # else :
    #     plt.scatter(houtx,houty,s=2,c=houtflux)
    #     plt.colorbar()
    #     plt.title('end of format bethermin')
    #     plt.show()

    print(maps[icol]['shead']['CRVAL1'], maps[icol]['shead']['CRVAL2'])
    #write everything to file for lenstool to ingest
    lensfile = (config.HOME + 'model/' + clusname + '/' + clusname + '_cat.cat')
    with open(lensfile,'w') as f :
        f.write('#REFERENCE 3 %.6f %.6f \n' %(maps[icol]['shead']['CRVAL1'], maps[icol]['shead']['CRVAL2']))
        for k in range(len(outmag)):
            f.write('%i %.3f %.3f 0.5 0.5 0.0 %0.6f %0.6f \n' \
                    %(k,houtx[k],houty[k],houtz[k],outmag[k]))
        f.close()

    truthtable = {'x': houtx, 'y': houty,
                  'z': houtz, 'mag': outmag}

    return retcat, truthtable

def get_gauss_beam(fwhm, pixscale, band, nfwhm=5.0, oversamp=1):
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
    clusname = 'a0370'
    resolution = 'nr'
    verbose = 1
    bolocam = None
    wave = [250.,350.,500.] # In units of um
    fwhm = [17.6, 23.9, 35.2]
    pixsize = [6.0, 8.33333, 12.0]
    fluxcut = 0
    import sys
    sys.path.append('../utilities')
    sys.path.append('../new_bethermin')
    sys.path.append('../source_handling')
    from clus_get_clusparams import clus_get_clusparams
    from clus_get_data import clus_get_data
    from genmap import genmap_gauss
    params, err = clus_get_clusparams(clusname,verbose=verbose)
    maps, err = clus_get_data(clusname=clusname, resolution=resolution, bolocam=bolocam, verbose=verbose)
    gm = genmap_gauss()
    sim_maps = gm.generate(0.25,verbose=True)
    clus_format_bethermin(0,sim_maps,maps,wave[0],clusname,pixsize[0],fluxcut=0,zzero=params['z'])
