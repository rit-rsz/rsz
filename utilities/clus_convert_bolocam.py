
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
import matplotlib.pyplot as plt
import math
from astropy.io import fits
import os
import pyfits

def clus_convert_bolocam(bolocam,verbose = 0,errmsg):

    bolocam.deconvolved_image = clus_dTtoDI(143,1e-6*bolocam.deconvolved_image)
    bolocam.deconvolved_image_smooth_trim_sn = \
        clus_dTtoDI(143,1e-6*bolocam.deconvolved_image_smooth_trim_sn)

    bolocamsize = (bolocam.deconvolved_image_noise_realizations).size
    for ib in range(bolocamsize):
#       This made use of reform but I am not sure that it is needed in python the way that it is used
        bolocam.deconvolved_image_noise_realizations[:,:,ib] = \
            clus_dTtoDI(143,1e-6*bolocam.deconvolved_image_noise_realizations)

reform(a,2,3)
	               arange(1,7).reshape(2,-1)
                    a.setshape(2,3

    return bolocam
