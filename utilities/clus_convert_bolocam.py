
################################################################################
# NAME : clus_convert_bolocam.py
# DATE STARTED : June 19, 2019
# AUTHORS : Dale Mercado
# PURPOSE : converts bolocam data in a way that it can be used in the SPIRE pipeline
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY : 3/6/202 - VLB updating to use new Bolocam template maps
################################################################################
import sys
from astropy.io import fits
sys.path.append('../sz')
from clus_dTtoDI import *
import numpy as np

def clus_convert_bolocam(bolocam, norm = None, clusname=None, verbose=0):

    if clusname == 'rxj1347' :
        # convert normalized map back to original dT
        bolocam_new = np.reshape([x*norm for x in bolocam.flatten()],(bolocam.shape[0],bolocam.shape[1]))
        bolocam_final = clus_dTtoDI(143,bolocam_new)

        return bolocam_final,None
    else :
        # Converts mK to MJy/sr
        bolocam[0]['deconvolved_image'] = clus_dTtoDI(143,1e-6*(bolocam[0]['deconvolved_image']))
        # Converts mK to MJy/sr
        bolocam[0]['deconvolved_image_smooth_trim_sn'] = \
            clus_dTtoDI(143,1e-6*bolocam[0]['deconvolved_image_smooth_trim_sn'])

        bolocamsize = bolocam[0]['deconvolved_image_noise_realizations'][0].shape

        # a 1-D array
        bolocam[0]['deconvolved_image_noise_realizations'] = \
            clus_dTtoDI(143,1e-6*bolocam[0]['deconvolved_image_noise_realizations'])
        bolocam[0]['deconvolved_image_noise_realizations'] = bolocam[0]['deconvolved_image_noise_realizations'].flatten()

        return bolocam,None
