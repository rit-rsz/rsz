
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

from math import *

def clus_dTtoDI(freq,dT):

    T_0     = 2.725                   #CMB temperature, K
    k_B     = 1.3806503e-23           #Boltzmann constant, J/K
    h       = 6.626068e-34            #Planck constant, J s
    HztoGHz = 1e9                     #Hz -> GHz
    m_e     = 5.11e2                  #electron mass keV
    c       = 2.99792458e8            #m/s

    nu = freq * HztoGHz
    x = (h * nu) / (k_B * T_0)
    fx = x * (exp(x) + 1) / (exp(x) - 1) - 4
    hx = x**4 * exp(x) / ((exp(x) - 1)**2)
    gx = hx * fx
    Iz = 2 * (k_B * T_0)**3 / (h * c)**2

    # Float converts to MJy as well as uK to K
    dI = (hx * (dT) / T_0)# * float(1e20)

    return dI
