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
sys.path.append('outputs')
from clus_get_clusparams import *
from clus_get_data import *
from clus_add_sziso_new import *
from clus_subtract_cat import *
from clus_subtract_xcomps import *
from clus_compute_rings import *
from clus_fitsz import *
from save_fitsz import *
from save_final_im import *
from bolocam_mask import *
import config
import argparse

class Catsrc():

    def __init__(self, clusname, isim=None, saveplot=1, maketf=0, sgen=None, verbose=1, resolution='nr', superplot=0, testflag=0, lense_only=0):
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
        clusters = ['a0370','a1689','a1835','a2218','a2219','a2390',
                  'cl0024','ms0451','ms1054','ms1358','rxj0152','rxj1347']
        self.beam = [get_spire_beam_fwhm('PSW'), #arcsec
                get_spire_beam_fwhm('PMW'),
                get_spire_beam_fwhm('PLW')]
        self.verbose = verbose
        self.saveplot = saveplot
        self.maketf = maketf
        self.sgen = sgen
        self.yin = config.yin[clusters.index(clusname)]
        self.tin = config.tin[clusters.index(clusname)]
        self.clusname = clusname
        self.nsim = int(isim) - 1
        self.resolution = resolution
        self.testflag = testflag
        self.superplot = superplot
        self.lense_only = lense_only
        self.dI = []
        self.data_retrieval()
        self.source_removal()
        self.data_analysis()
        # compile all figures and save to monolithic plot
        # save_final_im(self.sgen,self.nsim,self.clusname,testflag=self.testflag)

    def data_retrieval(self):
        """
        Purpose : the purpose of this function is to retrieve data from the
        files for the galaxies clusters and put them into a maps object
        Inputs : None
        Returns : None
        Class variables : Passes on self.maps (list of 3 dictionaries)
        """

        if self.verbose:
            print('Welcome to SZ fitter v 1.0 Python Version')
            print('Fetching cluster parameters')
        #fetch parameters for our cluster from a master csv file.
        params, err = clus_get_clusparams(self.clusname,self.nsim,verbose=self.verbose)
        if err:
            print(err)
            if self.verbose:
                print('clus_get_clusparams exited with error: ' + err)
            exit()

        if self.verbose:
            print('Fetching SPIRE maps')

        # fetch data from fits files and put them into a maps object.
        maps, err = clus_get_data(self.clusname,self.nsim,verbose=self.verbose,sgen=self.sgen,nsim=self.nsim, testflag=self.testflag)
        if err:
            if self.verbose:
                print('clus_get_data exited with error: ' + err)
            exit()

        # bolocam_mask(maps)

        # Add the sz effect into the simmulated clusters
        if self.sgen is not None and not self.lense_only:
            maps, err, dI = clus_add_sziso_new(maps,self.nsim,clusname=self.clusname,yin=self.yin,tin=self.tin,params=params,verbose=self.verbose, testflag=self.testflag,nsim=self.nsim,saveplot=self.saveplot)
            self.dI = dI
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

        maps, err = clus_subtract_cat(self.maps, self.dI, self.nsim, sgen=self.sgen, verbose=self.verbose, saveplot=self.saveplot, superplot=self.superplot)
        err = None

        # for i in range(3):
        #     self.maps[i]['srcrm'] = fits.getdata(config.OUTPUT + 'pcat_residuals/' + self.maps[i]['name'] + '_resid_' + self.maps[i]['band'] + '_' + str(self.nsim) + 'lense.fits')
        #     data_file = config.HOME + 'Lensing/lense_template_' + self.maps[i]['band'] + '.fits'
        #     lense_model = fits.getdata(data_file)
        #     self.maps[i]['srcrm'] = self.maps[i]['srcrm'] - lense_model

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
            print('Subtracting correlated components')

        maps, err = clus_subtract_xcomps(maps, sgen=self.sgen, verbose=self.verbose, superplot=self.superplot, saveplot=self.saveplot, nsim=self.nsim)
        if err:
            if self.verbose:
                print('clus_subtract_xcomps exited with error: ' + err)
            exit()

        if self.verbose:
            print('Saving processed images')

        self.maps = maps

    def data_analysis(self):
        """
        Purpose : this is the part of the code that actually does the data analysis for the pipeline
        Inputs : None
        Outputs : None
        Class Variables : self.maps gets passed. and passed out.
        """
        if self.verbose:
            print('Computing radial averages!')

        radave = clus_compute_rings(self.maps,self.params,30.0,sgen=self.sgen,verbose=self.verbose, superplot=self.superplot, saveplot=self.saveplot, nsim=self.nsim, testflag=self.testflag, lense_only=self.lense_only)

        if self.verbose:
            print('Computing Beta model fit.')
        #
        # if self.clusname == 'ms0451':
        #     maxlim = 300
        # else:
        #     maxlim = 450

        if self.lense_only == 0 :
            fit = clus_fitsz(radave, self.params,self.beam, sgen=self.sgen, superplot=self.superplot, saveplot=self.saveplot, nsim=self.nsim)
            fit = np.asarray(fit)
            increment = fit[:,0]
            offsets = fit[:,1]
            print('increment :',increment)
            print('offsets :', offsets)

            err = save_fitsz(increment, offsets, radave, self.params, sgen=self.sgen, verbose=self.verbose, nsim=self.nsim)
            if err:
                if self.verbose:
                    print('clus_save_szfits exited with error: ' + err)
                exit()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-run", help="This runs real or sim map through the pipeline for a single map on a specific cluster.", nargs=4, metavar=('clusname', 'sgen', 'nsim', 'resolution'))
    args = parser.parse_args()
    if args.run:
        print(args.run)
        clusname = args.run[0]
        sgen = args.run[1]
        nsim = args.run[2]
        resolution = args.run[3]
        catsrc = Catsrc(clusname, sgen=sgen, isim=nsim, resolution=resolution)
