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
# from clus_sz_template import * #used for transfer functions
sys.path.append('utilities')
from config import * #(this line will give me access to all directory variables)
from clus_get_tfs import *
from get_clusparams import *
from get_simmaps import *
sys.path.append('source_handling')
from subtract_cat import *
from subtract_xcomps import *
from get_data import *
from get_xid import *
sys.path.append('reduc')
from get_cats import *
sys.path.append('sz')
from add_sziso import *

class Catsrc():

    def __init__(self, clusname, saveplot=1, cattype="24um", savecat=0,
                 savemap=0, maketf=0, simmap=0, nsim=0, s2n=3, verbose=1, resolution='nr'):
        # possibly how we would get around defining these terms, not positive
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
        self.setup()

    def setup(self):
        if self.simmap > 0 and self.nsim == 0:
            if self.verbose:
                print('simmap set but nsim not supplied! Aborting')
            exit()

        if self.simmap == 0:
            self.nsim = None

        if self.verbose:
            print('Welcome to SZ fitter v 1.0 python')

        if self.saveplot:
            # Save the plots to an output ps file
            pass

        beam = [get_spire_beam_fwhm('PSW'), #arcsec
                get_spire_beam_fwhm('PMW'),
                get_spire_beam_fwhm('PLW')]

        if self.verbose:
            print('Fetching cluster parameters')
        params, err = get_clus_params(self.clusname,verbose=self.verbose)

        if err:
            if self.verbose:
                print('clus_get_clusparams exited with error: ' + err)
            exit()

        if self.verbose:
            print('Fetching SPIRE maps')

        # Get data now covers both simmulated and real data
        maps, err = get_data(self.clusname,verbose=self.verbose,simmap=self.simmap,nsim=self.nsim)
        if err:
            if self.verbose:
                print('clus_get_data exited with error: ' + err)
            exit()

        # Add the sz effect into the simmulated clusters
        if self.simmap:
            # First maps section is legacy
            # maps, err = get_simmaps(self.clusname,nsim=self.nsim, simflag=self.simmap, verbose=self.verbose)
            # if err:
            #     if self.verbose:
            #         print('clus_get_simmaps exited with error: ' + err)
            #     exit()
            maps, err = add_sziso(maps,yin=self.yin, tin=self.tin,verbose=self.verbose)
            if err:
                if self.verbose:
                    print('clus_add_sziso exited with error: '+ err)
                exit()

        ncols = len(maps)
        print(' ')

        if self.verbose:
            print('Fetching transfer functions')
        # ignore for now as this is only like a 2% correction and we are way off
        # if not self.maketf and not self.simmap:
        #     tf_maps, err = get_tfs(self.clusname)
        #     if err:
        #         tf_maps, err = get_data(self.clusname)
        #         ncols = len(tf_maps)
        #         for i in range(ncols):
        #             sztm = clus_sz_template(maps[i], params, verbose=self.verbose)
        #             tf_maps[i]['xclean'] = sztm #maybe need to change this
        #             tf_maps[i]['error'] = np.tile(1.0, tf_maps[i]['astr']['NAXIS'])
        #             tf_maps[i]['calfac'] = 1.0 / tf_maps[i]['JY2MJy']

        if self.verbose:
            print('Fetching regression catalogs')
        cat, err = get_cats(self.clusname,self.cattype,maps,savecat=self.savecat,
                                  savemap=self.savemap, simmap=self.simmap, nsim=self.nsim, s2n=self.s2n,
                                  verbose=self.verbose, resolution=self.resolution) #args need to be figured out
        if err:
            if self.verbose:
                print('clus_get_cats exited with error: ' + err)
            exit()

        if self.verbose:
            print('Band merging catalogs')

        xid, err = get_xid(maps, cat, savemap=self.savemap, simmap=self.simmap, verbose=self.verbose)
        if err:
            if self.verbose:
                print('clus_get_xid exited with error: ' + err)
            exit()

        if self.verbose:
            print('Regressing and subtracting catalogs')

        maps, err = subtract_cat(maps, xid, verbose=self.verbose)
        if err:
            if self.verbose:
                print('clus_subtract_cat exited with error: ' + err)
            exit()

        if self.verbose:
            print('Generating residual source mask')

        if not self.clusname:
            if self.verbose:
                print('Require a string array of cluster names as input, aborting!')

        #this is commented out because the function hasn't been made yet.
        #The idea is to use residual mask to do some manual masking, but we haven't
        #encountered the need to do that in our pipeline
        # residual_mask, err = clus_residual_mask(maps,verbose=self.verbose)
        # if err:
        #     if self.verbose:
        #         print('clus_residual_mask exited with error: ' + err)
        #         exit()

        if self.verbose:
            print('Subtracting correlated componenets')

        subtracted_comps, err = clus_subtract_xcomps(maps, simflag=self.simmap, verbose=self.verbose)
        if err:
            if self.verbose:
                print('clus_subtract_xcomps exited with error: ' + err)
            exit()

        if self.verbose:
            print('Saving processed images')

        err = clus_save_data(maps,yin=self.yin, tin=self.tin, simflag=self.simmap, verbose=self.verbose)
        if err:
            if self.verbose:
                print('clus_save_data exited with error: ' + err)
            exit()

        if self.verbose:
            print('Computing radial averages')

        radave, err = clus_compute_rings(maps,params,30.0,verbose=self.verbose)
        if err:
            if self.verbose:
                print('clus_compute_rings exited with error: ' + err)
        if self.simmap == None:
            tfave, err = clus_compute_rings(tf_maps, params, 30.0, verbose=self.verbose)
            if err:
                if self.verbose:

                    print('clus_compute_rings exited with error: ' + err)
                exit()

        # radave[2].fluxbin[0] = np.nan

        if self.verbose:
            print('Computing beta model fit.')

        if self.simmap:
            if self.clusname == 'ms0451':
                maxlim = 300
            else:
                maxlim = 450

            fit, err = clus_fitsz(radave,params,beam,offsets,maxlim=maxlim,minlim=0,superplot=superplot,
                                    verbose=self.verbose) #args need to be figued out when we write this function
            increment = fit[1,:] #don't know if this is the same as [1,*] in idl
            offsets = fit[0,:]
            if err:
                if self.verbose:
                    print('clus_fitsz exited with error: ' + err)
                exit()

            err = clus_save_szfits(increment, offsets, radave, params, simflag=self.simmap, verbose=self.verbose, outname='szout_' + str(nsim))

        if not self.simmap: #again not really sure if this is right.
            if self.clusname == 'ms0451':
                maxlim = 300
            else:
                maxlim = 450

            fit, err = clus_fitsz(radave,params,maxlim=maxlim,minlim=0,verbose=self.verbose) #args need to be worked out when we write the function
            increment = fit[1,:]
            offsets = fit[0,:]
            if err:
                if self.verbose:
                    print('clus_fitsz exited with error: ' + err)
                exit()

            err = clus_save_szfits(increment, offsets, radave, params, simflag=self.simmap, verbose=self.verbose)

        if err:
            if self.verbose:
                print('clus_save_szfits exited with error: ' + err)
            exit()

        if self.saveplots:
            # Save the plots to an output ps file
            pass


        return



if __name__ == '__main__':
    catsrc = Catsrc('a0370', verbose=1, cattype='PSW',simmap=2,nsim=200)
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
