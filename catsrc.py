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
import numpy as np
from math import *
import sys
sys.path.append('utilities')
sys.path.append('source_handling')
sys.path.append('reduc')
sys.path.append('sz')
from clus_get_clusparams import *
from clus_get_data import *
from clus_add_sziso import *
from clus_subtract_cat import *
from clus_subtract_xcomps import *
from clus_compute_rings import *
from clus_fitsz import *
from save_fitsz import *

class Catsrc():

    def __init__(self, clusname, saveplot=1, maketf=0, sgen=None, nsim=0, verbose=1, resolution='nr', superplot=1,
                 yin=False, tin=False, testflag=0):
        """
        initializing function for catsrc class
        Purpose: read in arguments to be passed to functions in catsrc.
        Inputs : Clusname   - the name of the galaxy cluster
                 saveplot   - flag to set if you want to save plots to a file
                 maketf     - make transfer function flag (currently not used)
                 sgen       - flag for looking at which rev of sims (None = real data, 0 = 0-99, 1 = 100 - 199, 2 = 200 - 299, etc.)
                 nsim       - the sim number
                 verbose    - verbosity flag
                 resolution - the resolution for our image
                 superplot  - flag to plot outputs to terminal during run
                 yin        - input compton y parameter
                 tin        - input temperature
                 testflag   - boolean flag for determining if you want to show tests for various parts of the pipeline
                              such as if the maps object is in the right order or if the rings
                              in compute rings are computed properly.
        Outputs :
                 szout - file containing the increment and offset as well as other information from
                 our analysis.
        """
        self.beam = [get_spire_beam_fwhm('PSW'), #arcsec
                get_spire_beam_fwhm('PMW'),
                get_spire_beam_fwhm('PLW')]
        self.verbose = verbose
        self.saveplot = saveplot
        self.maketf = maketf
        self.sgen = sgen
        self.yin = yin
        self.tin = tin
        self.clusname = clusname
        self.nsim = nsim
        self.resolution = resolution
        self.testflag = testflag
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
        if self.sgen is not None and self.nsim == 0:
            if self.verbose:
                print('sgen set but nsim not supplied! Aborting')
            exit()

        if self.sgen is None:
            self.nsim = 0

        if self.verbose:
            print('Welcome to SZ fitter v 1.0 Python Version')
            print('Fetching cluster parameters')
        #fetch parameters for our cluster from a master csv file.
        params, err = clus_get_clusparams(self.clusname,verbose=self.verbose)
        if err:
            if self.verbose:
                print('clus_get_clusparams exited with error: ' + err)
            exit()

        if self.verbose:
            print('Fetching SPIRE maps')

        # fetch data from fits files and put them into a maps object.
        maps, err = clus_get_data(self.clusname,verbose=self.verbose,sgen=self.sgen,nsim=self.nsim, testflag=self.testflag)

        if err:
            if self.verbose:
                print('clus_get_data exited with error: ' + err)
            exit()

        # Add the sz effect into the simmulated clusters
        if self.sgen is not None:
            maps, err = clus_add_sziso(maps,yin=self.yin, tin=self.tin,params=params,verbose=self.verbose, testflag=self.testflag)

        if err:
            if self.verbose:
                print('clus_add_sziso exited with error: '+ err)
            exit()

        #transfer function is not in use currently.
        # if self.verbose:
        #     print('Fetching transfer functions')
        # ignore for now as this is only like a 2% correction and we are way off
        # if not self.maketf and not self.sgen:
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

        maps, err = clus_subtract_cat(self.maps, verbose=self.verbose, saveplot=self.saveplot, nsim=self.nsim, superplot=self.superplot)
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

        maps, err = clus_subtract_xcomps(maps, sgen=self.sgen, verbose=self.verbose, superplot=self.superplot, saveplot=self.saveplot, nsim=self.nsim)
        if err:
            if self.verbose:
                print('clus_subtract_xcomps exited with error: ' + err)
            exit()

        if self.verbose:
            print('Saving processed images')

        #this is likely depreciated with the adittion of the saveplot flag.
        #it saves our maps object, but I don't think this is particularly useful.
        # err = clus_save_data(maps,yin=self.yin, tin=self.tin, sgen=self.sgen, verbose=self.verbose)
        # if err:
        #     if self.verbose:
        #         print('clus_save_data exited with error: ' + err)
        #     exit()


        self.maps = maps
        return maps

    def data_analysis(self):
        """
        Purpose : this is the part of the code that actually does the data analysis for the pipeline
        Inputs : None
        Outputs : None
        Class Variables : self.maps gets passed. and passed out.
        """
        if self.verbose:
            print('Computing radial averages!')

        radave = clus_compute_rings(self.maps,self.params,30.0,verbose=self.verbose, superplot=self.superplot, saveplot=self.saveplot, nsim=self.nsim)  #should superplot be a flag in catsrc?
        #unclear at the moment why we need to have two different calls to compute_rings
        # if self.sgen == None:  # don't see the difference between if sgen == 0 and if not sgen ??
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

        fit = clus_fitsz(radave, self.params,self.beam, superplot=self.superplot, saveplot=self.saveplot, nsim=self.nsim)
        fit = np.asarray(fit)
        increment = fit[:,0]
        print(increment)
        offsets = fit[:,1]


        #this is to do a specific case for ms0451, but I think we are cutting this cluster so I don't know if we need this.
        # if self.sgen is None:
        #     if self.clusname == 'ms0451':
        #         maxlim = 300
        #     else:
        #         maxlim = 450
        #
        #         fit = clus_fitsz(radave, params, self.beam)
        #         increment = fit[1,:]
        #         offsets = fit[0,:]

        err = save_fitsz(increment, offsets, radave, self.params, sgen=self.sgen, verbose=self.verbose, nsim=self.nsim)
        if err:
            if self.verbose:
                print('clus_save_szfits exited with error: ' + err)
            exit()
        self.offsets = offsets
        self.increment = increment
        return fit



if __name__ == '__main__':
    catsrc = Catsrc('a1689', verbose=1, sgen=2,nsim=200, superplot=0, tin=4, testflag=1)
