
################################################################################
# NAME : get_simmaps.py
# DATE STARTED : June 11, 2019
# AUTHORS : Victoria Butler & Dale Mercado
# PURPOSE : Much like get data this returns maps for the simmulated data
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
from astropy.io.ascii import read
from astropy.table import Column
# tab = read('spectrum.csv')
import csv
from collections import defaultdict


def get_simmaps():

    # Once I get this running I will put this into the get_simmaps inputs
    verbose = 1 if not verbose else verbose = verbose
    simflag = 1 if not simflag else simflag = simflag
    sb = 0 if not sb else sb = sb
    xc = 0 if not xs else xc = xc
