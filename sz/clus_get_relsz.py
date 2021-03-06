
################################################################################
# NAME : clus_get_relsz.py
# DATE STARTED : June 18, 2019
# AUTHORS : Victoria Butler & Dale Mercado
# PURPOSE : Given some nominal SZ inputs, this returns the expected SZ
#      effect signal.  Output is in MJy sr^-1.
# INPUTS :
#        nu  - the range of frequncies requested [fltarr, GHz]
#        y   - the Compton y parameter, usually O(10^-4) [float, unitless]
#        T_e - the electron temperature of the gas [float, keV]
#        vpec = cluster peculiar velocity, usually 0.0 [float, m/s]
#        ngrid = how finely sampled the frequency space should be [int]
#
# OUTPUTS : rSZ intensity at given SPIRE band frequency [MJy/sr]
# NOTES :
#
# REVISION HISTORY :
################################################################################
import scipy.io
import numpy as np
import sys,os
from subprocess import check_output
sys.path.append('../utilities')
from clus_get_lambdas import *
import config
from math import *
from scipy.interpolate import interp1d
from scipy import interpolate
sys.path.append('../lookup')
import matplotlib.pyplot as plt
##########################################################################################################
def clus_get_relsz(isim,nu,band,y=0,te=0,vpec=0.0, ngrid=100):

        # constants
        # '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
        # tau = []
        tau     = float(y) * (420.0 / float(te)**0.9) # optical depth
        T_0     = 2.725                   #CMB temperature, K
        k_B     = 1.3806503e-23           #Boltzmann constant, J/K
        h       = 6.626068e-34            #Planck constant, J s
        c       = 2.99792458e8            #m/s
        GHztoHz = 1.0e9                   #GHz -> Hz
        vpecp   = vpec * 1e3 / c            # peculiar velocity
        # '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
        # create input paramters for run_SZpack or read from existing params file
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
          '1e-6', # accuracy of numerical int
          ' ', # space
          '4',
          '2',
          ' ',
          '../lookup', # path for Output
          '.txt'] #filename extension

         # to accomodate multithreading with each band
        if band == 'BOLOCAM' :
            np.savetxt(config.HOME + 'lookup/szparams_%s.txt'%(band),sz_params,"%s",newline='\n')
            out = check_output([config.ROOT + 'SZpack.v1.1.1/run_SZpack','CNSNopt',config.HOME + 'lookup/szparams_%s.txt'%(band)])
        else :
            np.savetxt(config.HOME + 'lookup/szparams_%s_%s.txt'%(isim,band),sz_params,"%s",newline='\n')
            out = check_output([config.ROOT + 'SZpack.v1.1.1/run_SZpack','CNSNopt',config.HOME + 'lookup/szparams_%s_%s.txt'%(isim,band)])

        # run & parse 3 columns of stdout from run_SZpack into frequency and SZ signal arrays
        output = out.decode("utf-8") #convert from bytes to unicode
        # print(output)
        # print(output[:output.index('x=')-1]) # print header to screen
        try:
            output = output[output.index('x=')-1:] # remove header information before parsing
        except ValueError:
            pass

        output = output.split('\n')
        # output = large array with sub arrays of each line, where line = [x= val1, val2, val3]

        # remove all strings but data, & convert data to floats
        xout = []
        JofXout = []
        for line in range(len(output)-1):
            vals = output[line].strip(' x=').split(' ')        # plt.plot(xout,JofXout)

            vals = [float(i) for i in vals]
            xout.append(vals[0]) # frequency
            JofXout.append(vals[2]) # sz intensity [MJy/sr]
            # second column has kompaneets equation sz intensity [x^3*Dn(x)]
            # see Birkinshaw(1998) , pg.19

        ''' If reading from an archived data file '''
        # file = open('/home/butler/rsz/lookup/SZ_CNSN_basis.dat')
        # output = [x.strip('\n').split(' ') for x in file.readlines()]
        # old_xout = []
        # old_JofXout = []
        # for i in range(len(output)):
        #     old_xout.append(float(output[i][0]))
        #     old_JofXout.append(float(output[i][1]))
        ###################################################################
        # find the SZ intensity at a specific SPIRE frequency
        bigJ = interpolate.interp1d(xout,JofXout, kind='linear')
        if band == 'ALL':
            band_names = ['BOLOCAM', 'PSW', 'PMW', 'PLW']
            x_vals = [ (3e5 /  clus_get_lambdas(bandi, True))  * GHztoHz * h / (k_B * T_0)  for bandi in band_names]
            deltaI = [float(bigJ(x_valsi)) for x_valsi in x_vals]
            return deltaI, x_vals, JofXout, xout, None

        else:
            thisx = nu * GHztoHz * h / (k_B * T_0) # dimensionless frequency


            deltaI = float(bigJ(thisx))
            #This is plots for testing to see if this code is working.
            plt.plot(xout,JofXout)
            plt.scatter(thisx,deltaI,color='orange')
            plt.title('RUN_SZPACK SZE @ %0.3f GHZ : %0.4f [MJy/sr]' %(nu,deltaI))
            plt.xlabel('Dimensionless Frequency')
            plt.ylabel('$\Delta$I_0 [MJy/sr]')
            plt.savefig(config.HOME + 'lookup/sz_test_%s.png'%(band))
            plt.clf()
            # print('SZE @ %0.2f GHZ : %0.4f [MJy/sr] %0.6f %0.4f' %(nu,deltaI,y,te)) # print current return

            #remove old config file ...
            # if band == 'BOLOCAM' :
            #     if os.path.exists(config.HOME + 'lookup/szparams_%s.txt'%(band)):
            #         os.remove(config.HOME + 'lookup/szparams_%s.txt'%(band))
            # else :
            #     if os.path.exists(config.HOME + 'lookup/szparams_%s_%s.txt'%(isim,band)):
            #         os.remove(config.HOME + 'lookup/szparams_%s_%s.txt'%(isim,band))

            return deltaI, thisx, JofXout, xout, None



########################################################################################################################
if __name__ == '__main__':
    yin_coeff = [2.50,1.91,2.26,3.99,1.36,2.42,1.59,1.90,3.99]
    yin = [x*1e-4 for x in yin_coeff]
    tin = [7.2,10.1,7.7,9.8,4.5,8.6,7.8,5.5,10.9]
    print(yin[0],tin[0])
    clus_get_relsz(600.0,yin[0],tin[0])
    # sz_wrapper(argv[1],argv[2],argv[3],argv[4])
