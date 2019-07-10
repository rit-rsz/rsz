################################################################################
# NAME : clus_sz_template
# DATE STARTED : June 18, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE : to create a sz templates
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS : maps - a dictionary created by get_data
#          params - a dictionary created by get_params
#           verbose - 1 - if you want error messages
#                     0 - if you don't want error messages
#
#
# OUTPUTS : szmap -
# REVISION HISTORY :
################################################################################


from math import *
import numpy as np

def clus_sz_template(maps, params, verbose = 1):
    #shead is the fits header for the map

    racent = maps['shead']['CRVAL1']
    decent = maps['shead']['CRVAL2']
    rapix  = maps['shead']['CRPIX1']
    dcpix  = maps['shead']['CRPIX2']

    pixsize = maps['pixsize']

    dra = 3600.0 * (racent - params['fidrad']) # arcseconds
    ddc = -3600.0 * (decent- params['fidded']) #arcseconds
    dra = dra * cos(radians(decent)) / pixsize
    ddc = ddc / pixsize
    sizemap = maps['astr']['NAXIS1']
    midmapx = maps['astr']['CRPIX'][0] - 1 #float(sizemap[0] / 2 ) ???
    midmapy = maps['astr']['CRPIX'][1] - 1 #float(sizemap[1] / 2 ) ???
    print(sizemap)
    szmap = np.zeros([sizemap[0], sizemap[1]], dtype=float)

    rad_c = params['rc'] / pixsize #this is in arcsec
    beta = params['beta']
    norm = 1.0

    for i in range(sizemap[0]): #this is the same as the idl pipeline but this seems weird to me.
        for j in range(sizemap[1]): #same thing
            rad = sqrt((float(i) - midmapx - dra)**2 +
            (float(j) - midmapy - ddc)**2)
            szmap[i][j] = norm * (1.0 + (rad / rad_c)**2)**((1.0 - 3.0 * beta) / 2.0)

    return szmap
