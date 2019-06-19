
################################################################################
# NAME : catsrc.py
# DATE STARTED : June 11, 2019
# AUTHORS : Victoria Butler & Dale Mercado
# PURPOSE : This procedure is the first attempt at translating the old
#             pipeline from IDL into a working python counterpart.
#             The purpose is to act as a driver program that has the
#             capability to fit the SZ effect in the HerMES/HLS clusters.
#
#   Below is what Mike used to descrive the IDL scripts function
# ;;  This procedure is a first stab at a generalized driver program for
# ;;  fitting the SZ effect in the HerMES/HLS clusters.  It is designed so
# ;;  that one gives it a map and catalog as input and it fits for all the
# ;;  sources in the catalog from the map.  Then it takes those fits and
# ;;  subtracts or masks them from the map.  Finally, it takes radial averages
# ;;  and fits for the SZ effect.
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#           nsims (number of simulations/files to read in)
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import scipy.io
import numpy as np
from config import * #(this line will give me access to all directory variables)
import matplotlib.pyplot as plt
import math

import sys
sys.path.append('utilities')

import get_clusparams

# def catsrc(clusname,saveplots,cattype, savecat,savemap,maketf,simmap,nsim,s2n,yin,tin,verbose,success,errmsg):
class catsrc:
    # From when carsrc was a def
    # (clusname,\
    # saveplots=1,\
    # cattype=0,\
    # savecat=0,\
    # savemap=0,\
    # maketf=0,\
    # simmap=0,\
    # nsim=0,\
    # s2n=3,\
    # yin=0,\
    # tin=0,\
    # verbose=1):
    # ,SUCCESS=success,ERRMSG=errmsg):

    if not clusname:
        print('Require a string array of cluster names as input, aborting!')

    def __init__():
        self.maps[]
        # possibly how we would get around defining these terms, not positive
        self.verbose = 1 if not verbose else verbose
        self.cattype = '24 um' if not cattype else cattype
        self.savecat = 0 if not savecat else SAVECAT
        self.savemap = 0 if not savemap else SAVEMAP
        self.saveplot = 1 if not saveplot else saveplot
        self.maketf = 0 if not maketf else MAKETF
        self.sinmap = 0 if not SIMMAP else SIMMAP
        self.s2n = 3 if not s2n else S2N
        self.yin = 0  if not yin else YIN
        self.tin = 0 if not tin else tin



    def
#   Below assign values for these variables if they are not defined
#   I think I may have consolidated this into the def of catsrc
#   If not we can change back to these
        if simmap > 0 and not nsim:
            print('simmap set but nsim not supplied!   Aborting')

    #
    # if simmap == 0:
    #     nsim = #fill with nans??

    print(' ')
    if verbose:
        print('Welcome to SZ fitter v 1.0')
    print(' ')

    # if saveplots:
        #what ever the equivalent to nuplot would be
        # not sure if this is just setting a location or creating the file

    #This needs more work to be finished
    if verbose:
        print('Fetching cluster parameters')
    params = get_clusparams(clusname)

    if not cgp_success:
        # Need to make this a single string
        errmsg = ('clus_get_clusparams exited with error' + cgp_errm)
        if verbose:
            print(errmsg)

    if verbose:
        print('Fetching SPIRE maps')
    if not simmap:
        maps = get_data()
        if not cgd_success:
            errmsg = ('get_data exited with error' + cgd_errmsg)
            if verbose:
                print(errmsg)






if __name__ == '__main__':
    catsrc(clusname = 'rxj1347', verbose = 5)
        # SAVEPLOTS=saveplots,\
        # CATTYPE=cattype,\
        # SAVECAT=savecat,\
        # SAVEMAP=savemap,\
        # MAKETF=maketf,\
        # SIMMAP=simmap,\
        # NSIM=nsim,\
        # S2N=s2n,\
        # YIN=yin,\
        # TIN=tin,\
        # VERBOSE=verbose,SUCCESS=success,ERRMSG=errmsg)
