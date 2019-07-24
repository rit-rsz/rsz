
################################################################################
# NAME : clus_get_relsz.py
# DATE STARTED : June 18, 2019
# AUTHORS : Victoria Butler
# PURPOSE : Given some nominal SZ inputs, this returns the expected SZ
#      effect signal.  Output is in W m^-2 Hz^-1 sr^-1.
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#        nu  - the range of frequncies requested [fltarr, GHz]
#        y   - the Compton y parameter, usually O(10^-4) [float, unitless]
#        T_e - the electron temperature of the gas [float, keV
#
#
# OUTPUTS :
# NOTES :
#
# REVISION HISTORY :
################################################################################

# Currently this is in a rough state, while I am trying to interperate the idl script
import scipy.io
import numpy as np
import sys,os
from subprocess import check_output
sys.path.append('../utilities')
from config import * #(this line will give me access to all directory variables)
from math import *
from astropy.io import fits
from scipy.interpolate import interp1d
from scipy import interpolate

##########################################################################################################
def sz_wrapper(nu,y,te=10.0, vpec=0.0, ngrid=100):


        # tau = [x * (420.0 / te**0.9) for x in y]
        tau = y[0] * (420.0 / te**0.9)
        T_0     = 2.725                   #CMB temperature, K
        k_B     = 1.3806503e-23           #Boltzmann constant, J/K
        h       = 6.626068e-34            #Planck constant, J s
        c       = 2.99792458e8            #m/s
        HztoGHz = 1.0e9                     #Hz -> GHz
        vpecp = vpec * 1e3 / c

        if os.path.isfile('/home/mercado/rsz/lookup/szparams.txt') :
            pass
        else :
            sz_params = [
             '0.1', # xmin
             '30.0', # xmax
             ngrid, # x grid n_pts
             ' ', # space
             tau, # optical depth
             te, # electron temperature
             ' ', # space
             vpecp, # peculiar velocity in c
             '1.0', # direction cos of cluster v
             '0.0', # v_pec of observer wrt CMB
             '1.0', # direction of cos of observer
             ' ', # space
             '10', # temperature of order
             '2', # beta order
             '1e-4', # accuracy of numerical int
             ' ', # space
             '/home/butler/rsz/lookup', # path for Output
             '.txt'] #filename extension
            np.savetxt('/home/mercado/rsz/lookup/szparams.txt',sz_params,"%s",newline='\n')

        # child = subprocess.Popen(['/usr/local/bin/run_SZpack CNSN /home/butler/rsz/lookup/szparams.txt'],stdout=subprocess.PIPE)
        # out,err = child.communicate()
        out = check_output(['/usr/local/bin/run_SZpack','CNSN','/home/mercado/rsz/lookup/szparams.txt'])

        file = open('/home/mercado/rsz/lookup/SZ_CNSN_basis.dat')
        output = [x.strip('\n').split(' ') for x in file.readlines()]
        xout = []
        JofXout = []
        for i in range(len(output)):
            xout.append(float(output[i][0]))
            JofXout.append(float(output[i][1]))
        print(xout[-1],JofXout[-1])

        thisx = nu * HztoGHz * h / (k_B * T_0)
        print(thisx)
        I0 = 2 * (k_B * T_0)**3 / (h * c)**2

        # bigJ = INTERPOL(JofXout,xout,thisx)
        bigJ = interpolate.interp1d(xout,JofXout, kind='linear')
        act_bigJ = bigJ(thisx)
        deltaI = I0 * act_bigJ
        print(act_bigJ,deltaI)

        return deltaI, None



########################################################################################################################
if __name__ == '__main__':
    sz_wrapper(600.0,yin,10.0,0.0,100)
