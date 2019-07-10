
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
# REVISION HISTORY :
################################################################################
import scipy.io
import numpy as np
# from config import * #(this line will give me access to all directory variables)
from math import *
from astropy.io import fits
import os
import sys
sys.path.append('../sz')
from clus_dTtoDI import *

def clus_convert_bolocam(bolocam, verbose=0):
    # Converts mK to MJy/sr
    # bolocam['deconvolved_image'].setflags(write=1)
    print((bolocam[0]['deconvolved_image_smooth_trim_sn'])*1e-6)
    bolocam[0]['deconvolved_image'] = clus_dTtoDI(143,1e-6*(bolocam[0]['deconvolved_image']))
    bolocam[0]['deconvolved_image_smooth_trim_sn'] = \
        clus_dTtoDI(143,1e-6*bolocam[0]['deconvolved_image_smooth_trim_sn'])

    bolocamsize = (bolocam[0]['deconvolved_image_noise_realizations']).size
    for ib in range(bolocamsize-1):
#       This made use of reform but I am not sure that it is needed in python the way that it is used
        bolocam[0]['deconvolved_image_noise_realizations'][:,:,ib] = \
            clus_dTtoDI(143,1e-6*bolocam[0]['deconvolved_image_noise_realizations'])

        reform(a,2,3)
        arange(1,7).reshape(2,-1)
        a.setshape(2,3)
    print(bolocam)
    return bolocam,None
