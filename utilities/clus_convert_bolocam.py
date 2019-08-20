
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
import matplotlib.pyplot as plt

def clus_convert_bolocam(bolocam, verbose=0):
    # Converts mK to MJy/sr
    bolocam[0]['deconvolved_image'] = clus_dTtoDI(143,1e-6*(bolocam[0]['deconvolved_image']))

    bolocam[0]['deconvolved_image_smooth_trim_sn'] = \
        clus_dTtoDI(143,1e-6*bolocam[0]['deconvolved_image_smooth_trim_sn'])

    # This was used when the following conversion was in a for loop.
    # Due to issues with being not writable I'm trying to skip it and modify afterwards
    bolocamsize = bolocam[0]['deconvolved_image_noise_realizations'][0].shape

    # for ib in range(bolocamsize[2]):
    # This still needs to be reformed inorder to matcht the idl version. It should be
    # a 1-D array
    bolocam[0]['deconvolved_image_noise_realizations'] = \
        clus_dTtoDI(143,1e-6*bolocam[0]['deconvolved_image_noise_realizations'])
    bolocam[0]['deconvolved_image_noise_realizations'] = bolocam[0]['deconvolved_image_noise_realizations'].flatten()

    # print('bolocamsize = ', bolocam[0]['deconvolved_image_noise_realizations'][0])

    return bolocam,None
