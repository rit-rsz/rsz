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

def clus_format_bethermin(icol,sim_map,maps,band,clusname,pixsize,fwhm,\
                          fluxcut=0,zzero=0,superplot=0,savemaps=0):


    # trimming num sources for lenstool limit
    msrc = 50000 - 1

    # 3,4,5 are 250,350,500 truthtables
    cat = sim_map[icol+3]
    nsrc = len(cat['fluxdens']) # pulls len of truthtable
    refx = maps[0]['shead']['CRPIX1']
    refy = maps[0]['shead']['CRPIX2']

    # massage data into new arrays
    xpos = sim_map[icol+3]['x']
    ypos = sim_map[icol+3]['y']
    zpos = sim_map[icol+3]['z']
    outflux = sim_map[-1]['fluxdens'][:,icol]

    if superplot :
        plt.scatter(xpos,ypos,s=2,c=outflux)
        plt.colorbar()
        plt.title('Bethermin SIM (pre-format)')
        plt.show()

    outx = [pixsize * (x - refx) for x in xpos]
    outy = [pixsize * (y - refy) for y in ypos]
    outz = [float(np.ceil(10.0 * z)) / 10.0 for z in zpos]

    plt.scatter(outx,outy,s=2,c=outflux)
    plt.colorbar()
    plt.title('Bethermin SIM (pre-format)')
    plt.savefig('format_top_%s.png' %(band))
    plt.clf()
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
        pixscale = [6.0, 8.33333, 12.0]
        cmap = np.zeros((300,300),dtype=np.float32)

        if band == 'PSW':
            xf = np.floor(x)
            yf = np.floor(y)
            nx, ny = cmap.shape
            np.place(xf, xf > nx-1, nx-1)
            np.place(yf, yf > ny-1, ny-1)
            for cx, cy, cf in zip(xf, yf, flux):
                cmap[int(cy), int(cx)] += cf  # Note transpose

        # Other bands, with pixel scale adjustment
        else :
            if band == 'PMW':
                posrescale = pixscale[0] / pixscale[1]
            if band == 'PLW':
                posrescale = pixscale[0] / pixscale[2]

            xf = np.floor(posrescale * x)
            yf = np.floor(posrescale * y)
            nx, ny = cmap.shape
            np.place(xf, xf > nx-1, nx-1)
            np.place(yf, yf > ny-1, ny-1)
            for cx, cy, cf in zip(xf, yf, flux):
                cmap[int(cy), int(cx)] += cf  # Note transpose
            del posrescale

        beam = get_gauss_beam(fwhm,pixsize,band)
        sim_map = convolve(cmap, beam, boundary='wrap')

        # plt.imshow(sim_map,extent=(0,300,0,300),clim=[0.0,0.15],origin=0)
        # plt.colorbar()
        # plt.title('clus_format_bethermin: Non-lensed map')
        # plt.show()

        hdx = fits.PrimaryHDU(maps[icol]['signal'],maps[icol]['shead'])
        sz = fits.PrimaryHDU(sim_map,hdx.header)
        sz.writeto(config.SIMBOX + 'nonlensedmap_' + clusname + '_' + band + '.fits',overwrite=True)

    # magnitude instead of flux in Jy
    outmag = [-2.5 * np.log10(x) for x in houtflux]

    if superplot:
        plt.scatter(houtx,houty,s=2,c=houtflux)
        plt.colorbar()
        plt.title('end of format bethermin')
        plt.show()

    # write everything to file for lenstool to ingest
    lensfile = (config.HOME + 'model/' + clusname + '/' + clusname + '_cat.cat')
    new_k = []
    with open(lensfile,'w') as f :
        f.write('#REFERENCE 3 %.6f %.6f \n' %(maps[icol]['shead']['CRVAL1'], maps[icol]['shead']['CRVAL2']))
        for k in range(len(outmag)):
            new_k.append(k)
            f.write('%i %.3f %.3f 0.5 0.5 0.0 %0.6f %0.6f \n' \
                    %(k,houtx[k],houty[k],houtz[k],outmag[k]))
        f.close()

    np.save('k_%s.npy' %(band),new_k)

    return retcat

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
