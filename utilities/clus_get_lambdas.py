
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




def clus_get_lambdas(band):
    # wavelen was defined as lambda in IDL but python treats that an internal variable
    if band == 'PSW':
        wavelen = 250
    if band == 'PMW':
        wavelen = 350
    if band == 'PLW':
        wavelen = 500
    if band == 'BOLOCAM':
        wavelen = 2140

    return wavelen
