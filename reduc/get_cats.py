################################################################################
# NAME : get_cats.py
# DATE STARTED : June 21, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE : The purpose of this function is to fetch the catalog for a
#           specific type.
#
#
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import config
from astropy.io import fits
from math import *
import numpy as np
from astropy.wcs import WCS as wcs


def get_cats(clusname, cattype, maps, nsim, simmap=0, s2n=3, resoltuion='fr', verbose=1, savecat=0, savemap=0,):
    '''
    Purpose : to collect the catalogs given a specifc band
    Inputs : clusname - name of the cluster
             cattype - type of catalog
             maps - a dictionary created by get_data
             nsims - the number of sims
             simmap - 0 = no simmap
                      1 = simmap
            s2n -
            resolution - the resolution
            verbose - 0 = no error messages
                      1 = error messages
            savecat - 0 = don't save the catalog
                      1 = save the catalog
            savemap - 0 = don't save the map
                      1 = save the map
    Outputs: A dictionary containing a list of ras and decs, fluxes, err, clustername,
            instrumentm and s2n.
    '''
    err = False
    if cattype == 'PSW':
        if verbose:
            print('Requested %s catalog generation creating catalog' % (cattype))
        catalog, err = make_plw_src_cat(clusname, nsim, s2n=s2n,resolution=resolution, savecat=savecat, savemap=savemap, simmap=simmap, verbose=verbose) # this function needs to be written :)
        if err:
            if verbose:
                print(err)
            return None, err

    elif cattype == 'PLW':
        if verbose:
            print('Requested %s catalog generation creating catalog' % (cattype))
        catalog, err = make_plw_src_cat(clusname, nsim, resolution=resolution, simmap=simmap, s2n=s2n, verbose=verbose, savecat=savecat, savemap=savemap)
        if err:
            if verbose:
                print(err)
            return None, err

    elif cattype == '24um':
        if verbose:
            print('Requested %s catalog generation creating catalog' % (cattype))
        catalog, err = make_mips_src_cat(clusname, maps, s2n=s2n, savecat=savecat, savemap=savemap, verbose=verbose) # this also needs to be written :)
        if err:
            if verbose:
                print(err)
            return None, err
        print(cattype)

    elif cattype == 'MFLR':
        if verbose:
            print('Requested %s catalog generation creating catalog' % (cattype))
        catalog, err = make_mflr_src_cat(clusname, resolution=resolution, s2n=s2n, savecat=savecat, savemap=savemap, verbose=verbose)
        if err:
            if verbose:
                print(err)
            return None, err

    else:
        err = 'cattype not recognized %s' % cattype
        if verbose:
            print(err)
        return None, err

    return catalog, err

def make_plw_src_cat(clusname, resolution, nsim, simmap=0, s2n=3, verbose=1, savecat=0, savemap=0):
    '''
    Purpose: to create a catalog with the PLW band.\
    Inputs: clusname - clustername
            resolution - resoltuion
            nsim - the number of sims.
            simmap - 0 = no simmap
                     1 = simmap
            s2n -
            verbose - 0 = no error messages
                      1 = error messages
            savecat - 0 = don't save the catalog
                      1 = save the catalog
            savemap - 0 = don't save the map
                      1 = save the map
    Outputs: cat - the catalog dictionary.
    '''
    err = False
    if s2n < 3 and verbose:
        print('WARNING: S/N requested less than 3, extraction may be suspect')

    if not simmap:
        min_corr = 0.5
        maps, err = get_data(clusname, resolution=resolution, verbose=verbose)
        if err:
            if verbose:
                print('clus_get_data exited with error:' + err)
            return None, err

    else:
        min_corr = 0.01
        maps, err = get_data(clusname, simflag=simflag, nsim=nsim, verbose=verbose)
        if err:
            if verbose:
                print('clus get data exited with error: ' + err)
            return None, err
    for i in range(len(maps)):
        if 'PLW' in maps[i]['band']:
            index = i
    dataPLW = maps[index]['signal']
    errPLW = maps[index]['error']
    headPLW = maps[index]['shead']
    psfPLW = maps[index]['psf']

    if verbose:
        print('Constructing PLW catalog')

    #starfinder here we are going to ignore that for now

    if savemap:
        if verbose:
            print('saving catalog debug maps')
        writefits(config.CLUSSBOX + 'make_PLW_src_cat_model.fits', data=modelPSW, header_dict=headPSW)
        data = np.empty(dataPLW.shape)
        for i in range(dataPLW.shape[0]):
            for j in range(dataPLW.shape[1]):
                data[i,j] = dataPLW[i,j] - modelPLW[i,j]
        writefits(config.CLUSSBOX + 'make_PLW_src_cat_mdiff.fits', data=data, header_dict=headPSW)
    astr = extast(headPLW)
        #more starfinder stuff computing RA and dec from x, y.

        #another call here to more stuff from starfinder
    ra_dec = wcs.wcs_pix2world(xPLW, yPLW, origin, ra_dec_order=True)
        #xPLW = a list of x coordinates
        #yPLW = a list of y coordinates
        #origin = Idk what to put for the origin
        #ra_dec is going to be a list of ra/dec pairs.

        # whpl = WHERE(fPLW/sigf GE s2n,count)
        # a = a[whpl]
        # d = d[whpl]
        # fPLW = fPLW[whpl]
        # sigf = sigf[whpl]

    if verbose:
        print('CUT S/N >= %s sources, kept %s stars' % (s2n, count))

    if savecat:
        catfile = config.CLUSDATA + 'catalogs/' + clusname + '_PLW.dat'
        if verbose:
            print('Saving catalog data to %s' % (catfile))
        file = open(catfile, 'w')
        file.write('Catalog extracted by make_PLW_src_cat')
        file.write('Created: ' + str(datetime.datetime.now()))
        file.write('Extracted from ' + imgfile)
        file.write('Extracted S/N >= ' + str(s2n))
        file.write('RA    DEC    FLUX    ERROR')
        for i in range(len(a)):
            myline = '9B' + str(a[i]).translate(str.maketrans('', '', string.whitespace)) + '9B' + \
                            str(d[i]).translate(str.maketrans('', '', string.whitespace)) + '9B' + \
                            str(fPSW[i]).translate(str.maketrans('', '', string.whitespace)) + '9B' + \
                            str(sigf[i]).translate(str.maketrans('', '', string.whitespace))
            file.write(myline)
        file.close()
    cat = {'ra': a,
           'dec' : d,
           'flux' : fPLW,
           'err' : sigf,
           'cluster' : clusname,
           'instrument': 'SPIRE PLW',
           's2n' : s2n}
    return cat, err

def make_mflr_src_cat(clusname, resolution='fr', s2n=3, savecat=0, savemap=0, verbose=1):
    '''
    Purpose: to create a catalog with the MFLR cattype.
    Inputs: clusname - clustername
            resolution - resoltuion
            nsim - the number of sims.
            simmap - 0 = no simmap
                     1 = simmap
            s2n -
            verbose - 0 = no error messages
                      1 = error messages
            savecat - 0 = don't save the catalog
                      1 = save the catalog
            savemap - 0 = don't save the map
                      1 = save the map
    Outputs: cat - the catalog dictionary.
    '''
    err = False
    if s2n < 2 and verbose:
        print('WARNING: S/N requested less than 2, extraction may be suspect')
    min_corr = 0.1
    maps, err = get_data(clusname, resolution=resolution, verbose=verbose)
    if err:
        if verbose:
            print('clus_get_data exited with error: ' + err)
        return None, err

    for i in range(len(maps)):
        if 'PSW' in maps[i]['band']:
            index = i
    newmaps = maps[i]

    filtmaps = clus_matched_filter(newmaps) #this function also needs to be written :) :)

    dataPSW = filtmaps['signal']
    errPSW = filtmaps['error']
    headPSW = filtmaps['shead']
    psfPSW = filtmaps['psw']

    if savemap:
        if verbose:
            print('Saving matched filter map')
            writefits(config.CLUSSBOX + 'make_MFLR_src_cat_filt.fits', data=dataPSW, header_dict=headPSW)

    if verbose:
        print('Constructing MFLR catalog')

    #ANOTHER call to starfinder, which we don't have so RIP

    if savemap:
        if verbose:
            print('Saving catalog debug maps.')
        writefits(config.CLUSSBOX + 'make_PSW_src_cat_model.fits', data=modelPSW, header_dict=headPSW)
        data = np.empty(dataPSW.shape)
        for i in range(dataPSW.shape[0]):
            for j in range(dataPSW.shape[1]):
                data[i,j] = dataPSW[i,j] - modelPSW[i,j]
        writefits(config.CLUSSBOX + 'make_PSW_src_cat_mdiff.fits', data=data, header_dict=headPSW)

    astr = extast(map)

    ra_dec = wcs.wcs_pix2world(xPSW, yPSW, origin, ra_dec_order=True)
    #xPLW = a list of x coordinates
    #yPLW = a list of y coordinates
    #origin = Idk what to put for the origin
    #ra_dec is going to be a list of ra/dec pairs.

    # whpl = WHERE(fPSW/sigf GE s2n,count)
    # a = a[whpl]
    # d = d[whpl]
    # fPSW = fPSW[whpl]
    # sigf = sigf[whpl]

    if verbose:
        print('Cut S/N >= %s sources, kept %s stars' % (s2n, count))

    if savecat:
        catfile = config.CLUSDATA + 'catalogs/' + clusname + '_MFLR.dat'
        if verbose:
            print('Saving catalog data to %s' % (catfile))
        file = open(catfile, 'w')
        file.write('Catalog extracted by make_MFLR_src_cat')
        file.write('Created :' + str(datetime.datetime.now()))
        file.write('Extracted from ' + imgfile)
        file.write('Extracted S/N >= ' + s2n)
        file.write('RA    DEC    FLUX    ERROR')
        for i in range(len(a)):
            myline = '9B' + str(a[i]).translate(str.maketrans('', '', string.whitespace)) + '9B' + \
                      str(d[i]).translate(str.maketrans('', '', string.whitespace)) + '9B' + \
                      str(fPSW[i]).translate(str.maketrans('', '', string.whitespace)) + '9B' + \
                      str(sigf[i]).translate(str.maketrans('', '', string.whitespace))
            file.write(myline)
        file.close()

    cat = {'ra': a,
           'dec' : d,
           'flux' : fPLW,
           'err' : sigf,
           'cluster' : clusname,
           'instrument': 'SPIRE PLW',
           's2n' : s2n}
    return cat, err

def make_psw_src_cat(clusname, resolution, nsim, s2n=3, savecat=0, savemap=0, simmap=0, verbose=1):
    '''
    Purpose: to create a catalog with the PSW band.\
    Inputs: clusname - clustername
            resolution - resoltuion
            nsim - the number of sims.
            simmap - 0 = no simmap
                     1 = simmap
            s2n -
            verbose - 0 = no error messages
                      1 = error messages
            savecat - 0 = don't save the catalog
                      1 = save the catalog
            savemap - 0 = don't save the map
                      1 = save the map
    Outputs: cat - the catalog dictionary.
    '''
    err = False
    if s2n < 3 and verbose:
        print('WARNING: S/N requested less than 3, extraction may be suspect')

    if not simmap:
        min_corr = 0.5
        maps, err = get_data(clusname, resolution=resolution, verbose=verbose)
        if err:
            if verbose:
                print('Clus_get_data exited with error: ' + err)
            return None, err
    else:
        min_corr = 0.1
        maps, err = get_simmaps(clusname, simflag=simflag, nsim=nsim,verbose=verbose)
        if err:
            if verbose:
                print('clus_get_data exited with error: ' + err)
            return None, err

    for i in range(len(maps)):
        if 'PSW' in maps[i]['band']:
            index = i

    dataPSW = maps[index]['signal']
    errPSW = maps[index]['error']
    headPSW = maps[index]['shead']
    psfPSW = maps[index]['psf']

    if verbose:
        print('constructing PSW catalog')

    #STARINFDER call to starfinder.

    if savemap:
        if verbose:
            print('Saving catalog debug maps')
        writefits(config.CLUSSBOX + 'make_PSW_src_cat_model.fits', data=modelPSW, header_dict=headPSW)
        data = np.empty(dataPSW.shape)
        for i in range(dataPSW.shape[0]):
            for j in range(dataPSW.shape[1]):
                data[i,j] = dataPSW[i,j] - modelPSW[i,j]
        writefits(config.CLUSSBOX + 'make_PSW_src_cat_mdiff.fits', data=data, header_dict=headPSW)

    astr = extast(headPSW)

    ra_dec = wcs.wcs_pix2world(xPSW, yPSW, origin, ra_dec_order=True)
    #xPSW = a list of x coordinates
    #yPSW = a list of y coordinates
    #origin = Idk what to put for the origin
    #ra_dec is going to be a list of ra/dec pairs.

    # whpl = WHERE(fPSW/sigf GE s2n,count)
    # a = a[whpl]
    # d = d[whpl]
    # fPSW = fPSW[whpl]
    # sigf = sigf[whpl]

    if verbose:
        print('Cut S/N >= %s sources, kept %s stars' % (s2n, count))

    if savecat:
        catfile = config.CLUSDATA + 'catalogs/' + clusname + '_PSW.dat'
        if verbose:
            print('Saving catalog data to %s .' % (catfile))
        file = open(catfile, 'w')
        file.write('Catalog extracted by make_PSW_src_cat')
        file.write('Created :' + str(datetime.datetime.now()))
        file.write('Extracted from ' + imgfile)
        file.write('Extracted S/N >= ' + str(s2n))
        file.write('RA    DEC    FLUX    ERROR')
        for i in range(len(a)):
            myline = '9B' + str(a[i]).translate(str.maketrans('', '', string.whitespace)) + '9B' + \
                            str(d[i]).translate(str.maketrans('', '', string.whitespace)) + '9B' + \
                            str(fPSW[i]).translate(str.maketrans('', '', string.whitespace)) + '9B' + \
                            str(sigf[i]).translate(str.maketrans('', '', string.whitespace))
            file.write(myline)
        file.close()

    cat = {'ra': a,
           'dec' : d,
           'flux' : fPLW,
           'err' : sigf,
           'cluster' : clusname,
           'instrument': 'SPIRE PLW',
           's2n' : s2n}
    return cat, err

def make_mips_src_cat(clusname, maps, s2n=3, savecat=0, savemap=0, verbose=1):
    '''
    Purpose: to create a catalog with the 24um cattype.\
    Inputs: clusname - clustername
            resolution - resoltuion
            nsim - the number of sims.
            simmap - 0 = no simmap
                     1 = simmap
            s2n -
            verbose - 0 = no error messages
                      1 = error messages
            savecat - 0 = don't save the catalog
                      1 = save the catalog
            savemap - 0 = don't save the map
                      1 = save the map
    Outputs: cat - the catalog dictionary.
    '''
    err = False
    if s2n < 3:
        print('WARNING: S/N requested less tahn 3, extraction may be suspect')

    imgfile = config.CLUSDATA + 'mips/' + clusname + '_24um_sig.fits'
    uncfile = config.CLUSDATA + 'mips/' + clusname + '_24um_unc.fits'
    psffile = config.CLUSDATA + 'mips/mips_psf/mips_24_3000K.fits'


    if imgfile:
        hdul = fits.open(imgfile)
        imgh2 = hdul
        imgh = hdul[0].header
        img = hdul[0].data
    else:
        err = 'Cannot find %s' %(imgfile)
        if verbose:
            print('Cannot find %s' % s (imgfile))
        return None, err
    if uncfile:
        hdul = fits.open(uncfile)
        unch = hdul[0].header
        unc = hdul[0].data
        #there are some operations here but i don't know what they do
    else:
        err = 'Cannot find %s' %(uncfile)
        if verbose:
            print(err)
        return None, err
    if psffile:
        hdul = fits.open(psffile)
        psfh = hdul[0].header
        psf = hdul[0].data
    else:
        err = 'Cannot find %s' %(psffile)
        if verbose:
            print(err)
        return None, err

    imgpix = sqrt(imgh['PXSCAL1'] * imgh['PXSCAL1'] +  imgh['PXSCAL2'] * imgh['PXSCAL2']) / sqrt(2.0)
    psfpix = psfh['PIXSCALE']

    #call to change image scale don't know what that does so :)

    #call to starfinder

    if savemap:
        if verbose:
            print('Saving catalog debug maps.')
        writefits(config.CLUSSBOX + 'make_mips_src_cat_model.fits', data=imgmodel, header_dict=imgh)
        data = np.empty(imgmodel.shape)
        for i in range(imgmodel.shape[0]):
            for j in range(imgmodel.shape[1]):
                data[i,j] = img[i,j] - imgmodel[i,j]
        writefits(config.CLUSSBOX + 'make_mips_src_cat_mdiff.fits', data=data, header_dict=imgh)

    astr = extast(imgh2[0])
    count = 0

    ra_dec = wcs.wcs_pix2world(x, y, origin, ra_dec_order=True)
    #x = a list of x coordinates
    #y = a list of y coordinates
    #origin = Idk what to put for the origin
    #ra_dec is going to be a list of ra/dec pairs.

    # whpl = WHERE(f/sf GE s2n,count)    # whpl = WHERE(f/sf GE s2n,count)
    # a = a[whpl]
    # d = d[whpl]
    # f = f[whpl]
    # sf = sf[whpl]
    # a = a[whpl]
    # d = d[whpl]
    # f = f[whpl]
    # sf = sf[whpl]
    if verbose:
        print('Cut S/N >= %s sources, kept %s stars' %(s2n, count))
    if maps:
        if verbose:
            print('Generating mask cuts')
        ncol = len(maps)
        for i in range(ncol):
            #call to HASTROM i don't think we have an equivalent lol
            # whpl = np.where(np.isfinite(unc_align) == False)
            # whpl2 = np.where(unc_align >= 0.2)
            # maps[i]['mask'][whpl] =  1
            # maps[i]['mask'][whpl2] = 1
            #call to write to a fits image here
            pass
    if savecat:
        catfile = config.CLUSDATA + 'catalogs/' + clusname + '_24um.dat'
        if verbose:
            print('Saving catalog data to %s' % (catfile))
        file = open(catfile, 'w')
        file.write('Catalog extracted by make_mips_src_cat')
        file.write('Created: ' + str(datetime.datetime.now()))
        file.write('Extracted from ' + imgfile)
        file.write('Extracted S/N >= ' + str(s2n))
        file.write('RA    DEC    FLUX    ERROR')
        for i in range(len(a)):
            myline = '9B' + str(a[i]).translate(str.maketrans('', '', string.whitespace)) + '9B' + \
                            str(d[i]).translate(str.maketrans('', '', string.whitespace)) + '9B' + \
                            str(f[i]).translate(str.maketrans('', '', string.whitespace)) + '9B' + \
                            str(sf[i]).translate(str.maketrans('', '', string.whitespace))
            file.write(myline)
        file.close()

    cat = {'ra': a,
           'dec' : d,
           'flux' : fPLW,
           'err' : sigf,
           'cluster' : clusname,
           'instrument': 'SPIRE PLW',
           's2n' : s2n}
    return cat, err


def extast(map):
    '''
    Purpose : gather information from a FITS header file.
    Inputs : map - a dictionary of information generated by get_data.py.
    Outputs : astr - a dictionary containing information from the header file.
    '''
    astr = {}
    try:
        cd11 = map.header['CD1_1']
        cd12 = map.header['CD1_2']
        cd21 = map.header['CD2_1']
        cd22 = map.header['CD2_2']
        pv1_1 = map.header['PV1_0']
        pv1_2 = map.header['PV1_1']
        pv1_3 = map.header['PV1_2']
        pv1_4 = map.header['PV1_3']
        pv1_5 = map.header['PV1_4']
    except KeyError:
        pass # i don't like the way this is coded probably have to change it later
    for keys in map.header.keys():
        if 'NAXIS' in keys:
            astr.update({keys : map.shape})
        if 'CD1_1' in keys:
            x =  np.array([[cd11, cd12], [cd21, cd22]])
            astr.update({'CD' : x})
        if 'CDELT' in keys:
            astr.update({keys : map.header[keys]})
        if 'CRPIX1' in keys:
            x = np.array([map.header['CRPIX1'], map.header['CRPIX2']])
            astr.update({'CRPIX' : x})
        if 'CTYPE1' in keys:
            x = np.array([map.header['CTYPE1'], map.header['CTYPE2']])
            astr.update({'CTYPE' : x})
        if 'CRVAL1' in keys:
            x = np.array([map.header['CRVAL1'], map.header['CRVAL2']])
            astr.update({'CRVAL' : x})
        if 'LONGPOLE' in keys:
            astr.update({keys : map.header[keys]})
        if 'LATPOLE' in keys:
            astr.update({keys : map.header[keys]})
        if 'PV1_0' in keys:
            x = np.array([pv1_1, pv1_2, pv1_3, pv1_4, pv1_5])
            astr.update({keys : x})
        return astr
