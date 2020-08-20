
################################################################################
# NAME : add_lambdas.py
# DATE STARTED : June 20, 2019
# AUTHORS : Dale Mercado
# PURPOSE : Adding in the SZ isomap
# EXPLANATION : This script assigns the wavelength in um to the given band.
# CALLING SEQUENCE :
# INPUTS : band defined as one of the four listed in this file.
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################




def clus_get_lambdas(band, center=False):
        #band_w = [2140, 262.0563,361.63 ,514]
    # wavelen was defined as lambda in IDL but python treats that an internal variable
    if center:
        if band == 'PSW':
            wavelen = 262
        elif band == 'PMW':
            wavelen = 361.63
        elif band == 'PLW':
            wavelen = 514
        elif band == 'BOLOCAM':
            wavelen = 2140
    else:
        if band == 'PSW':
            wavelen = 250
        elif band == 'PMW':
            wavelen = 350
        elif band == 'PLW':
            wavelen = 500
        elif band == 'BOLOCAM':
            wavelen = 2140

    return wavelen
