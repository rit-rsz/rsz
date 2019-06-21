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

def get_cats(clusname, cattype, maps,resolution, simmap, nsim, s2n, verbose=1, savecat=0, savemap=0,):
    if cattype == 'PSW':
        if verbose:
            print('Requested %s catalog generation creating catalog' % (cattype))
        catalog, err = make_plw_src_cat(args) # this function needs to be written :)
        if err:
            if verbose:
                print(err)
            return None, err

    if cattype == 'PLW':
        if verbose:
            print('Requested %s catalog generation creating catalog' % (cattype))
        catalog, err = make_plw_src_cat(args)
        if err:
            if verbose:
                print(err)
            return None, err

    if cattype == '24um':
        if verbose:
            print('Requested %s catalog generation creating catalog' % (cattype))
        catalog, err = make_mips_src_cat(args) # this also needs to be written :)
        if err:
            if verbose:
                print(err)
            return None, err

    if cattype == 'MFLR':
        if verbose:
            print('Requested %s catalog generation creating catalog' % (cattype))
        catalog, err = make_mflr_src_cat(args)
        if err:
            if verbose:
                print(err)
            return None, err

    else:
        err = 'cattype not recognized'
        if verbose:
            print(err)
        return None, err

    return catalog

def make_plw_src_cat(clusname, resolution, nsim, simmap=0, s2n=3, verbose=1, savecat=0, savemap=0):
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
                print('clus get data exited with error: ' err)
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
        #writefits is called here need to make a replacement for this type of function
        #don't know what dictionary to add to the fits file.
    astr = extast(headPLW)
        #more starfinder stuff computing RA and dec from x, y.

        #another call here to more stuff from starfinder
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
        #code that then writes to a file

    cat = {'ra': a,
           'dec' : d,
           'flux' : fPLW,
           'err' : sigf,
           'cluster' : clusname,
           'instrument': 'SPIRE PLW',
           's2n' : s2n}
    return cat

def make_mflr_src_cat(clusname, resolution='fr', s2n=3, savecat=0, savemap=0, verbose=1):
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
            print('SAving matched filter map')
        #code to save stuff to a fits file still need to figure this out.

    if verbose:
        print('Constructing MFLR catalog')

    #ANOTHER call to starfinder, which we don't have so RIP

    if savemap:
        if verbose:
            print('Saving catalog debug maps.')
        #another call to write fits so RIP

    astr = extast(map)

    #converting form xy to celestial coordinate system don't have that code yet
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
        #code to then save the data to a file

    cat = {'ra': a,
           'dec' : d,
           'flux' : fPLW,
           'err' : sigf,
           'cluster' : clusname,
           'instrument': 'SPIRE PLW',
           's2n' : s2n}
    return cat

def make_psw_src_cat(clusname, resolution, nsim, s2n=3, savecat=0, savemap=0, simmap=0, verbose=1):
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
        #code to write to a fits file

    astr = extast(headPSW)

    #converting from xy to celestial coordinates again.
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
        #code to save to a file

    cat = {'ra': a,
           'dec' : d,
           'flux' : fPLW,
           'err' : sigf,
           'cluster' : clusname,
           'instrument': 'SPIRE PLW',
           's2n' : s2n}
    return cat

def make_mips_src_cat(clusname, maps, s2n=3, savecat=0, savemap=0, verbose=1):
    if s2n < 3:
        print('WARNING: S/N requested less tahn 3, extraction may be suspect')

    imgfile = config.CLUSDATA + 'mips/' + clusname + '_24um_sig.fits'
    uncfile = config.CLUSDATA + 'mips/' + clusname + '_24um_unc.fits'
    psffile = config.CLUSDATA + 'mips/mips_psf/mips_24_3000k.fits'

    if imgfile:
        hdul = fits.open(imgfile)
        imgh2 = hdul
        imgh = hdul[1].header
        img = hdul[1].data
    else:
        err = 'Cannot find %s' %(imgfile)
        if verbose:
            print('Cannot find %s' %s (imgfile))
        return None, err
    if uncfile:
        hdul = fits.open(uncfile)
        unch = hdul[1].header
        unc = hdul[1].data
        #there are some operations here but i don't know what they do
    else:
        err = 'Cannot find %s' %(uncfile)
        if verbose:
            print(err)
        return None, err
    if psffile:
        hdul = fits.open(psffile)
        psfh = hdul[1].header
        psf = hdul[1].data
    else:
        err = 'Cannot find %s' %(psffile)
        if verbose:
            print(err)
        return None, err

    imgpix = math.sqrt(imgh['PXSCAL1'] * imgh['PXSCAL1'] +  imgh['PXCAL2'] * imgh['PXCAL2']) / math.sqrt(2.0)
    psfpix = psfh['PIXSCALE']

    #call to change image scale don't know what that does so :)

    #call to starfinder

    if savemap:
        if verbose:
            print('Saving catalog debug maps.')
        #code to write to a fits file

    astr = extast(imgh2)

    #call to convert from xy to celestial coordinate system
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
            whpl = np.where(np.isfinite(unc_align) == False)
            whpl2 = np.where(unc_align >= 0.2)
            maps[i]['mask'][whpl] =  1
            maps[i]['mask'][whpl2] = 1
            #call to write to a fits image here
    if savecat:
        catfile = config.CLUSDATA + 'catalogs/' + clusname + '_24um.dat'
        if verbose:
            print('Saving catalog data to %s' % (catfile))
        #code to save data to a file

    cat = {'ra': a,
           'dec' : d,
           'flux' : fPLW,
           'err' : sigf,
           'cluster' : clusname,
           'instrument': 'SPIRE PLW',
           's2n' : s2n}
    return cat


def extast(map):
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
