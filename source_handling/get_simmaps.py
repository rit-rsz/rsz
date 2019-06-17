
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

    # Need to figure out where thismap in the idl script is being grabbed from in the idl script
    # Once I get this running I will put this into the get_simmaps inputs
    verbose = 1 if not verbose else verbose = verbose
    simflag = 1 if not simflag else simflag = simflag
    sb = 0 if not sb else sb = sb
    xc = 0 if not xs else xc = xc

    if simflag > 0 and not nsim:
        errmsg = 'simmap set but nsim not supplied! Aborting'
        if verbose:
            return errmsg
            # Can't use a goto like in idl or other older languages

    if simflag = 0:
        nsim = np.nan
        cols = ['PSW','PMW','PLW']
        ncols = len(cols)
        # How do want to set maps =
        # In idl this is where maps is then made into a pointer

    if simflag == 1:
        for i in range(len(col)-1):
            mapfile = CLUSDATA + 'sim_clusters/' + clusname + '_' + cols[i] + '.sav'
            mapfilename = scipy.io.readsav(mapfile, python_dict = True)
            # !!!!! This tries to use a calfac that is included from somewhere.
            # thismap.calfac = 1/(thismap.calfac * (thismap.pixsize)**2)
            # This is a different calibration factor than the one we have worked with before
            if sb:
                overmap = CLUSDATA + 'sim_cluster/' + clusname +\
                          '_' + cols[i] + '_sb.fits'
                overmap = fits.open(overmap)
                # Is this something that I have to define in the innit file???
                thismap.xclean = overmap
                whpl =
