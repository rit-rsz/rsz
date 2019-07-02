################################################################################
# NAME : matched_filter.py
# DATE STARTED : July 2, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE : wrapper script for applying matched_filter.pro
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import numpy as np

def matched_filter():
    conf_def = [5.8, 6.3, 6.8] * 1*10**3 # Jy / Beam
    fwhm_nom = [18.0, 25.0, 36.0]

    map = mapin['signal']
    noise = mapin['error']

    #sanitize maps -- set missing pix to largeval * max(noise)
    badind = []
    goodind = []
    for i in range(noise.shape[0]):
        for j in range(noise.shape[1]):
            if np.isfinite(noise[i,j]) == 0 or noise[i,j] == 0:
                badind.append([i,j])
            else:
                goodind.append([i,j])
    if len(badind) > 0:
        largeval = 100.0
        map[badind] = 0.0
        noise[badind] = largeval * np.amax(noise[goodind])

    #find band index
    if mapin['band'] == 'PSW':
        bandind = 0
    elif mapin['band'] == 'PMW':
        bandind = 1
    elif mapin['band'] == 'PLW':
        bandind = 2

    #set fwhm
    fwhm_arc = fwhm_nom[bandind]

    #fwhm in pixels
    fwhm_pix = fwhm_arc / mapin['pixsize']

    #MATCHED_FILTER(ARGS)

def matched_filter_pro():
    
