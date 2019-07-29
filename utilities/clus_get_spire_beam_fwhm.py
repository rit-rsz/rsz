################################################################################
# NAME : get_spire_beam_fwhm.py
# DATE STARTED : June 18, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE : Fetches the full-width half-maximum for a given band.
# EXPLANATION : Units are arcseconds
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################

def get_spire_beam_fwhm(band):
    if band == 'PSW':
        beamFWHM = 18.0
        return beamFWHM
    elif band == 'PMW':
        beamFWHM = 25.0
        return beamFWHM
    elif band == 'PLW':
        beamFWHM = 36.0
        return beamFWHM
    elif band == 'BOLOCAM':
        beamFWHM = 144.0
        return beamFWHM
    else:
        print('Unknown band:', band)
