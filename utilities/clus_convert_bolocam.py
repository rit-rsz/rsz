
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
from math import *
from astropy.io import fits
import os
import sys
sys.path.append('../sz')
from clus_dTtoDI import *
import matplotlib.pyplot as plt

def clus_convert_bolocam(bolocam, verbose=0):
    # Converts mK to MJy/sr
    bolocam[0]['deconvolved_image'] = clus_dTtoDI(143,1e-6*(bolocam[0]['deconvolved_image']))
    # Converts mK to MJy/sr
    bolocam[0]['deconvolved_image_smooth_trim_sn'] = \
        clus_dTtoDI(143,1e-6*bolocam[0]['deconvolved_image_smooth_trim_sn'])

    bolocamsize = bolocam[0]['deconvolved_image_noise_realizations'][0].shape

    bolocam[0]['deconvolved_image_noise_realizations'] = \
        clus_dTtoDI(143,1e-6*bolocam[0]['deconvolved_image_noise_realizations'])
    bolocam[0]['deconvolved_image_noise_realizations'] = bolocam[0]['deconvolved_image_noise_realizations'].flatten()


    return bolocam,None
