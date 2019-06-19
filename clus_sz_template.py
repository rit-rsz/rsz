import math as m


def clus_sz_template(maps, params, verbose = 1):
    #maps should be a fits header array object
    racent = maps['CRVAL1'] #i think this is the way to do it, but not sure.
    decent = maps['CRVAL2']
    rapix  = maps['CRPIX1']
    dcpix  = maps['CRPIX2']
    pixsize = maps.pixsize # don't know if this one will work.

    dra = 3600.0 * (racent - params.fidrad) # arcseconds
    ddc = -3600.0 * (decent- params.fidded) #arcseconds
    dra = dra * m.cos(m.radians(decent)) / pixsize
    ddc = ddc / pixsize

    sizemap = maps.astr.naxis #again don't know if this will work.
    sizemapx = maps.astr.crpix[0] - 1 #float(sizemap[0] / 2 ) ???
    sizemapy = maps.astr.crpix[1] - 1 #float(sizemap[1] / 2 ) ???

    szmap = np.zeroes(sizemap, dtype=float)

    rad_c = params.rc / pixsize #this is in arcsec
    beta = params.beta
    norm = 1.0

    for i in range(sizemap[0]): #not sure if len(sizemap[0]) or just sizemap[0]
        for j in range(sizemap[1]): #same thing
            rad = m.sqrt((float(i) - midmapx - dra)**2 +
            (float(j) - midmapy - ddc)**2)
            szmap[i][j] = norm * (1.0 + (rad / rad_c)**2)**((1.0 - 3.0 * beta) / 2.0)

    return szmap
