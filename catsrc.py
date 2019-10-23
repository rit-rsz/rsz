################################################################################
# NAME : catsrc.py
# DATE STARTED : June 11, 2019
# AUTHORS : Victoria Butler & Dale Mercado & Benjamin Vaughan
# PURPOSE : This procedure is the first attempt at translating the old
#             pipeline from IDL into a working python counterpart.
#             The purpose is to act as a driver program that has the
#             capability to fit the SZ effect in the HerMES/HLS clusters.
#
#   Below is what Mike used to describe the IDL scripts function
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
import matplotlib.pyplot as plt
from math import *
import sys
# from clus_sz_template import *
sys.path.append('utilities')
from config import * #(this line will give me access to all directory variables)
# from clus_get_tfs import * #this isn't being used right now
from clus_get_clusparams import *
# from clus_pcat_setup import *
# This is no longer being used
# from get_simmaps import *
sys.path.append('source_handling')
from clus_subtract_cat import *
from clus_subtract_xcomps import *
from clus_get_data import *
import config
from clus_get_xid import *
sys.path.append('reduc')
from clus_get_cats import *
sys.path.append('sz')
from clus_add_sziso import *
from clus_compute_rings import *
from clus_fitsz import *
from save_fitsz import *
sys.path.append('multiband_pcat')
# from multiband_pcat import *

class Catsrc():

    def __init__(self, clusname, saveplot=1, cattype="24um", savecat=0,
                 savemap=0, maketf=0, simmap=0, nsim=0, s2n=3, verbose=1, resolution='nr', superplot=1):
        """
        initializing function for catsrc class
        Purpose: read in arguments to be passed to functions in catsrc.
        Inputs : Clusname   - the name of the galaxy cluster
                 saveplot   - flag to set if you want to save plots to a file
                 cattype    - the catalog type (this might be changed in the future when pcat is integrated)
                 savecat    - flag to save the catalogs (this might change with pcat)
                 savemap    - flag to save the map object
                 maketf     - make transfer function flag (currently not used)
                 simmap     - flag for looking at which rev of sims
                 nsim       - the sim number
                 s2n        - s2n the signal to noise ratio
                 verbose    - verbosity flag
                 resolution - the resolution for our image
        Outputs : None
        """
        self.beam = [get_spire_beam_fwhm('PSW'), #arcsec
                get_spire_beam_fwhm('PMW'),
                get_spire_beam_fwhm('PLW')]
        self.verbose = verbose
        self.cattype = cattype
        self.savecat = savecat
        self.savemap = savemap
        self.saveplot = saveplot
        self.maketf = maketf
        self.simmap = simmap
        self.s2n = s2n
        self.yin = yin
        self.tin = tin
        self.clusname = clusname
        self.nsim = nsim
        self.resolution = resolution
        self.superplot = superplot
        self.data_retrieval()
        self.source_removal()
        self.data_analysis()

    def data_retrieval(self):
        """
        Purpose : the purpose of this function is to retrieve data from the
        files for the galaxies clusters and put them into a maps object
        Inputs : None
        Returns : None
        Class variables : Passes on self.maps (list of 3 dictionaries)
        """
        if self.simmap > 0 and self.nsim == 0:
            if self.verbose:
                print('Simmap set but nsim not supplied! Aborting')
            exit()

        if self.simmap == 0:
            self.nsim = np.nan

        if self.verbose:
            print('Welcome to SZ fitter v 1.0 Python Version')

        if self.saveplot:
            #this was used to plot some data and graph, but is depreciated in the python
            #version maybe we want to add it back?
            pass



        if self.verbose:
            print('Fetching cluster parameters')
        params, err = clus_get_clusparams(self.clusname,verbose=self.verbose)
        if err:
            if self.verbose:
                print('clus_get_clusparams exited with error: ' + err)
            exit()

        if self.verbose:
            print('Fetching SPIRE maps')

        # This step now is for both sims and real data
        maps, err = clus_get_data(self.clusname,verbose=self.verbose,simmap=self.simmap,nsim=self.nsim)

        if err:
            if self.verbose:
                print('clus_get_data exited with error: ' + err)
            exit()

        # Add the sz effect into the simmulated clusters
        if self.simmap:
            maps, err = clus_add_sziso(maps,yin=self.yin, tin=self.tin,params=params,verbose=self.verbose)
            # plt.imshow(maps[2]['signal'])
            # plt.show()

        if err:
            if self.verbose:
                print('clus_add_sziso exited with error: '+ err)
            exit()

        #transfer function is not in use currently.
        ncols = len(maps)
        # if self.verbose:
        #     print('Fetching transfer functions')
        # ignore for now as this is only like a 2% correction and we are way off
        # if not self.maketf and not self.simmap:
        #     tf_maps, err = get_tfs(self.clusname)
        #     if err:
        #         tf_maps, err = clus_get_data(self.clusname)
        #         ncols = len(tf_maps)
        #         for i in range(ncols):
        #             sztm = clus_sz_template(maps[i], params, verbose=self.verbose)
        #             tf_maps[i]['xclean'] = sztm #maybe need to change this
        #             tf_maps[i]['error'] = np.tile(1.0, tf_maps[i]['astr']['NAXIS'])
        #             tf_maps[i]['calfac'] = 1.0 / tf_maps[i]['JY2MJy']

        self.maps = maps
        self.params = params
        return




    def source_removal(self):
        """
        Purpose : the purpose of this function is to act as a wrapper for the source removal script that
        we are using. That used to be XID but most of this code has been depreciated.
        Inputs : self.maps get passed by reference
        Outputs : None
        Class variables : self.maps gets passed in and set on the way out
        """

        if self.verbose:
            print('Regressing and subtracting catalogs')

        maps, err = clus_subtract_cat(self.maps, verbose=self.verbose)
        if err:
            if self.verbose:
                print('clus_subtract_cat exited with error: ' + err)
            exit()

        if self.verbose:
            print('Generating residual source mask')

        if not self.clusname:
            if self.verbose:
                print('Require a string array of cluster names as input, aborting!')

        #This is commented out because the function hasn't been made yet.
        #The idea is to use residual mask to do some manual masking, but we haven't
        #encountered the need to do that in our pipeline
        # residual_mask, err = clus_residual_mask(maps,verbose=self.verbose)
        # if err:
        #     if self.verbose:
        #         print('clus_residual_mask exited with error: ' + err)
        #         exit()

        if self.verbose:
            print('Subtracting correlated componenets')

        maps, err = clus_subtract_xcomps(maps, simflag=self.simmap, verbose=self.verbose, superplot=self.superplot)
        if err:
            if self.verbose:
                print('clus_subtract_xcomps exited with error: ' + err)
            exit()

        if self.verbose:
            print('Saving processed images')

        # err = clus_save_data(maps,yin=self.yin, tin=self.tin, simflag=self.simmap, verbose=self.verbose)
        # if err:
        #     if self.verbose:
        #         print('clus_save_data exited with error: ' + err)
        #     exit()

        if self.simmap == 0:
            if self.verbose:
                print('Computing radial averages nsim=200')
        self.maps = maps
        return maps

    def data_analysis(self):
        """
        Purpose : this is the part of the code that actually does the data analysis for the pipeline
        Inputs : None
        Outputs : None
        Class Variables : self.maps gets passed. and passed out.
        """
        radave = clus_compute_rings(self.maps,self.params,30.0,verbose=self.verbose, superplot=self.superplot)  #should superplot be a flag in catsrc?
        #unclear at the moment why we need to have two different calls to compute_rings
        # if self.simmap == None:  # don't see the difference between if simmap == 0 and if not simmap ??
        #     tfave, err = clus_compute_rings(tf_maps, params, 30.0, verbose=self.verbose)
        #     if err:
        #         if self.verbose:
        #             print('clus_compute_rings exited with error: ' + err)
        #         exit()

        # radave[2].fluxbin[0] = np.nan #i guess this is right??

        if self.verbose:
            print('Computing Beta model fit.')

        if self.clusname == 'ms0451':
            maxlim = 300
        else:
            maxlim = 450

        fit = clus_fitsz(radave, self.params,self.beam, superplot=self.superplot) #args need to be figued out when we write this function
        fit = np.asarray(fit)
        increment = fit[:,0]
        print(increment)
        offsets = fit[:,1]

        if not self.simmap: #again not really sure if this is right.
            if self.clusname == 'ms0451':
                maxlim = 300
            else:
                maxlim = 450

                fit = clus_fitsz(radave, params, self.beam)
                increment = fit[1,:]
                offsets = fit[0,:]

            # increment = increment / tfamp # i don't think tfamp is defined?
        err = save_fitsz(increment, offsets, radave, self.params, simflag=self.simmap, verbose=self.verbose, outname='szout_')
        if err:
            if self.verbose:
                print('clus_save_szfits exited with error: ' + err)
            exit()
        self.offsets = offsets
        self.increment = increment
        # if self.saveplots:
        #     pass
            #we want to save some plots from the pipeline here.
        return fit



if __name__ == '__main__':
    catsrc = Catsrc('a0370', verbose=1, cattype='PSW',simmap=2,nsim=200, superplot=0)
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
