
################################################################################
# NAME : clus_dttodi.py
# DATE STARTED : June 19, 2019
# AUTHORS : Dale Mercado
# PURPOSE : converts dT to dI
# EXPLANATION : This is needed to convert the bolocam data into the same
#               format as the SPIRE data.
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import scipy.io
import numpy as np
import math
# from config import * #(this line will give me access to all directory variables)

def clus_dttodi(freq,dT):

    T_0     = 2.725                   #CMB temperature, K
    k_B     = 1.3806503e-23           #Boltzmann constant, J/K
    h       = 6.626068e-34            #Planck constant, J s
    HztoGHz = 1e9                     #Hz -> GHz
    m_e     = 5.11e2                  #electron mass keV
    c       = 2.99792458e8            #m/s

    nu = freq * HztoGHz

    x = (h * nu) / (k_b * T_0)

    fx = x * (exp(x) + 1.) / (exp(x) - 1.) - 4.d0
    hx = x^4 * exp(x) / ((exp(x) - 1.)^2)
    gx = hx * fx
    Iz = 2.d0 * (k_B * T_0)^3 / (h * c)^2

    dI = Iz * hx * dT / T_0
#   Below is a remnent for what I assume is for a wider data set, left in incase this needs to be turned on
#   if freq < 217:
#       dI = (-1)*dI

    return dI
