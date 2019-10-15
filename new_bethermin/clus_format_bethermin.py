
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
import sys
sys.path.append('../utilities')
import config
import numpy as np


def clus_format_bethermin(icol,sim_map,maps,band,clusname,pixsize,\
                          fluxcut=0,zzero=0):


    # initialize variables

    msrc = 50000 - 1

    # 3,4,5 are 250,350,500 truthtables
    cat = sim_map[icol+3]
    nsrc = len(cat['fluxdens']) # pulls len of truthtable
    refx = maps[icol]['shead']['CRPIX1']
    refy = maps[icol]['shead']['CRPIX2']
    print(refx,refy)

    # massage data into new arrays
    xpos = sim_map[icol+3]['x']
    ypos = sim_map[icol+3]['y']
    zpos = sim_map[icol+3]['z']
    print(xpos[0:10],ypos[0:10])
    outx = [pixsize * (x - refx) for x in xpos]
    outy = [pixsize * (y - refy) for y in ypos]
    # outx = [pixsize * x for x in xpos]
    # outy = [pixsize * y for y in ypos]
    print(outx[0:10],outy[0:10])
    outz = [float(np.ceil(10.0 * z)) / 10.0 for z in zpos]
    outflux = sim_map[-1]['fluxdens'][:,icol]

    # lets do some fluxcuts
    if fluxcut > 0 :
        print('Cutting input catalog at flux %s' %(fluxcut))
        for i in range(len(outflux)):
            if outflux[i] <= fluxcut :
                np.delete(outflux,i)
                np.delete(outx,i)
                np.delete(outy,i)
                np.delete(outz,i)
        nsrc = len(outflux)

    print(len(outflux))
    savex = []
    savey = []
    savez = []
    savef = []
    if zzero > 0 :
        for j in range(len(outflux)):
            if outz[j] <= zzero :
                savef.append(outflux[j])
                np.delete(outflux,j)
                savex.append(outx[j])
                np.delete(outx,j)
                savey.append(outy[j])
                np.delete(outy,j)
                savez.append(outz[j])
                np.delete(outz,j)
        retcat = {'x':savex,'y':savey,'z':savez,'f':savef}
        nsrc = len(outflux)

    print(len(outflux))
    # sort according to brightness due to lenstool limitations
    outx = [x for _,x in sorted(zip(outflux,outx), reverse=True)]
    outy = [y for _,y in sorted(zip(outflux,outy), reverse=True)]
    outz = [z for _,z in sorted(zip(outflux,outz), reverse=True)]
    outflux = sorted(outflux, reverse=True)

    # truncate to the msrc brightest sources
    if msrc < nsrc :
        outflux = outflux[0:msrc]
        outx = outx[0:msrc]
        outy = outy[0:msrc]
        outz = outz[0:msrc]
    print(len(outflux))
    # now sort according to z
    outflux = [f for _,f in sorted(zip(outz,outflux))]
    outx = [x for _,x in sorted(zip(outz,outx))]
    outy = [y for _,y in sorted(zip(outz,outy))]
    outz = sorted(outz)

    # magnitude instead of flux in Jy
    outmag = [-2.5 * np.log10(x) for x in outflux]

    # write everything to file for lenstool to ingest
    lensfile = (config.HOME + 'model/' + clusname + '/' + clusname + '_cat.cat')
    with open(lensfile,'w') as f :
        f.write('#REFERENCE 3 %.6f %.6f \n' %(maps[icol]['shead']['CRVAL1'], maps[icol]['shead']['CRVAL2']))
        for k in range(len(outmag)):
            f.write('%i %.3f %.3f 0.5 0.5 0.0 %0.6f %0.6f \n' \
                    %(k,outx[k],outy[k],outz[k],outmag[k]))
        f.close()

    return retcat

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
