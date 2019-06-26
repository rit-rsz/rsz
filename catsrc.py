
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
from config import * #(this line will give me access to all directory variables)
import matplotlib.pyplot as plt
from math import *
import sys
from clus_sz_template import *
sys.path.append('utilities')
from clus_get_tfs import *
from get_clusparams import *
sys.path.append('source_handling')
from get_data import *
import config
from get_xid import *
sys.path.append('reduc')
from get_cats import *
# def catsrc(clusname,saveplots,cattype, savecat,savemap,maketf,simmap,nsim,s2n,yin,tin,verbose,success,errmsg):
class Catsrc():

    def __init__(self, clusname,nsim=0, verbose=1, cattype="24um", savecat=0, savemap=0, saveplot=1,
                 maketf=0, simmap=0, s2n=3,yin=0,tin=0):
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
        self.setup()

    def setup(self):
        if self.simmap > 0 and len(self.nsim) == 0:
            if self.verbose:
                print('simmap set but nsim not supplied! Aborting')
            exit()

        if self.simmap == 0:
            self.nsim = np.nan #maybe this isn't what is being set.

        if self.verbose:
            print('Welcome to SZ fitter v 1.0')

        if self.saveplot:
            #call to nuplot idk what to put here instead
            pass

        ringw = 18.0 #arcseconds
        calfac = (pi/ 180.0) * (1.0 / 3600)*(1.0 / 3600) * (pi / 4 * log(2.0)) * 1*10**6
        PMWthres = 5*10**-3
        PLWthres = 8*10**-3

        ncols = 3.0

        #beam = [get_spire_beam_fwhm('PSW'), get_spire_beam_fwhm('PMW'), get_spire_beam_fwhm('PLW')]

        if self.verbose:
            print('Fetching cluster parameters')
        params, err = get_clus_params(self.clusname,verbose=self.verbose)
        if err:
            if self.verbose:
                print('clus_get_clusparams exited with error: ' + err)
            exit()

        if self.verbose:
            print('Fetching SPIRE maps')

        if not self.simmap:
            maps, err = get_data(self.clusname,verbose=self.verbose)
            if err:
                if self.verbose:
                    print('clus_get_data exited with error: ' + err)
                exit()
        else:
            maps, err = clus_get_simmaps(self.clusname, simflag=self.simmap, verbose=self.verbose, nsim=self.nsim)
            if err:
                if self.verbose:
                    print('clus_get_simmaps exited with error: ' + err)
                exit()
            maps, err = clus_add_sziso(maps,yin=self.yin, tin=self.tin,verbose=self.verbose)
            if err:
                if self.verbose:
                    print('clus_add_sziso exited with error: '+ err)
                exit()
        ncols = len(maps)

        if self.verbose:
            print('Fetching transfer functions')
        if not self.maketf and not self.simmap:
            tf_maps, err = get_tfs(self.clusname)
            if err:
                tf_maps, err = get_data(self.clusname)
                ncols = len(tf_maps)
                for i in range(ncols):
                    sztm = clus_sz_template(maps[i], params, verbose=self.verbose)
                    tf_maps[i]['xclean'] = sztm #maybe need to change this
                    tf_maps[i]['error'] = np.tile(1.0, tf_maps[i]['astr']['NAXIS'])
                    tf_maps[i]['calfac'] = 1.0 / tf_maps[i]['JY2MJy']
        if self.maketf:
            if self.verbose:
                print('Making SZ template map.')
            sztm, err = clus_sz_template(maps[0], params, verbose=self.verbose)
            if err:
                if self.verbose:
                    print('clus_sz_template exited with error: '+ err)
                exit()

            sztemplatemapfile = config.CLUSBOS + self.lusname + '_sztemplate.fits'
            #this part creates a fits header file from data and i need to properly work out how to do this

            intfmaps = np.empty(ncols)
            for i in range(ncols):
                intfmaps[i] = maps[i].file[0]

            conffiles = config.SMAP_PATH + 'map_making/createmap/conffiles/' + clusname + '_tf.conf' #config.SMAP_PATH needs to be filled in
            itest = smap_transfun_pipline_v4(args) #args need to be worked out when writing this function

            tf_maps, err = clus_get_tfs(self.clusname)
            if err:
                if self.verbose:
                    print('Cannot find created TF files, Aborting.')
                exit()

        if self.verbose:
            print('Fetching regression catalogs')
        cats, err = get_cats(self.clusname,self.cattype,maps,savecat=self.savecat,
                                  savemap=self.savemap, simmap=self.simmap, nsim=self.nsim, s2n=self.s2n,
                                  verbose=self.verbose) #args need to be figured out
        if err:
            if self.verbose:
                print('clus_get_cats exited with error: ' + err)
            exit()

        if self.verbose:
            print('Band merging catalogs')
        xid, err = get_xid(maps, cats, savemap=self.savemap, simmap=self.simmap, verbose=self.verbose)
        if err:
            if self.verbose:
                print('clus_get_xid exited with error: ' + err)
            exit()

        if self.verbose:
            print('Regressing and subtracting catalogs')

        maps, err = clus_subtract_cat(cats, maps, xid, verbose=self.verbose)
        if err:
            if self.verbose:
                print('clus_subtract_cat exited with error: ' + err)
            exit()

        if self.verbose:
            print('Generating residual source mask')

        if not self.clusname:
            if self.verbose:
                print('Require a string array of cluster names as input, aborting!')
        residual_mask, err = clus_residual_mask(maps,verbose=self.verbose)
        if err:
            if self.verbose:
                print('clus_residual_mask exited with error: ' + err)
                exit()

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

        if self.simmap == 0:
            if self.verbose:
                print('Computing radial averages')

            radave, err = clus_compute_rings(maps,params,30.0,verbose=self.verbose)
            if err:
                if self.verbose:
                    print('clus_compute_rings exited with error: ' + err)
        if self.simmap == None:  # don't see the difference between if simmap == 0 and if not simmap ??
            tfave, err = clusnames_compute_rings(tf_maps, params, 30.0, verbose=self.verbose)
            if err:
                if self.verbose:    # From when carsrc was a def

                    print('clus_compute_rings exited with error: ' + err)
                exit()

        radave[2].fluxbin[0] = np.nan #i guess this is right??

        if self.verbose:
            print('Computing beta model fit.')

        if self.clusname == 'ms0451':
            maxlim = 300
        else:
            maxlim = 450

        fit, err = clus_fitsz(args) #args need to be figued out when we write this function
        increment = fit[1,:] #don't know if this is the same as [1,*] in idl
        print('***** Increment value in catsrc ******', increment)
        offsets = fit[0,:]
        if err:
            if self.verbose:
                print('clus_fitsz exited with error: ' + err)
            exit()

        if not self.simmap: #again not really sure if this is right.
            if self.clusname == 'ms0451':
                maxlim = 300
            else:
                maxlim = 450
            print(maxlim)

            fit, err = clus_fitsz(args) #args need to be worked out when we write the function
            increment = fit[1,:]
            print('***** Increment Value 2 in catsrc ******', increment)
            offsets = fit[0,:]
            if err:
                if self.verbose:
                    print('clus_fitsz exited with error: ' + err)
                exit()
            increment = increment / tfamp # i don't think tfamp is defined?

            err = clus_save_szfits(increment, offsets, radave, params, simflag=self.simmap, verbose=self.verbose)

        else:
            err = clus_save_szfits(increment, offsets, radave, params, simflag=self.simmap, verbose=self.verbose, outname='szout_' + str(nsim))
        if err:
            if self.verbose:
                print('clus_save_szfits exited with error: ' + err)
            exit()

        if self.saveplots:
            pass
            #do something i don't know what UNPLOT does

        return #idk if we need to return anything here lol.



if __name__ == '__main__':
    catsrc = Catsrc('rxj1347', verbose=1)
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