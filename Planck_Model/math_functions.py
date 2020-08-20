################################################################################
# NAME : make_map.py
# DATE STARTED : June 25, 2020
# AUTHORS : Benjamin Vaughan
# PURPOSE : This set of code contains the mathematical functions for calculating
# cirrus emission using a modified blackbody from https://arxiv.org/pdf/1312.1300.pdf
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import numpy as np
from math import *

def blackbody_func(T, nu):
    '''
    This is a standard Planck Law
    Inputs : T (float array) - Temperatures of dust
             Nu (float) - the frequency of interest
    Outputs : bb (float array) - an array of results from Planck Law at given T and Nu
    '''
    k       = 1.3806503e-23           #Boltzmann constant, J/K
    h       = 6.626068e-34            #Planck constant, J s
    c       = 2.99792458e8            #m/s
    const =  2 * h * nu**3 / c**2
    exponent = np.divide(h * nu / k, T)

    denom = np.exp(exponent) - 1
    bb = np.divide(const, denom)
    return bb

def calc_intensity(beta, tau,T, nu):
    '''
    This takes the Planck Law function and does the calculation for the MBB
    Inputs : beta (float array) - beta parameter for power law of nu dependence at a given pixel
             tau (float array) - the optical depth at a given pixel
             T (float array) - the temperature at a given pixel
             Nu (float) - the frequency of interest
    Outputs: I_nu (float array) - Intensity of each pixel in SI units
    '''
    nu_const = nu / 353e9 #GHz
    power = np.asarray([nu_const**bi for bi in beta])
    bb = blackbody_func(T, nu)
    I = np.multiply(tau, bb)
    I_nu = np.multiply(I, power)
    return I_nu
