################################################################################
# NAME : clus_sim_background.py
# DATE STARTED : September 27, 2019
# AUTHORS : Dale Mercado & Victoria Butler
# PURPOSE : The wrapper for Conley's Bethermin simmulation generator.
#           Fluxcut is in Jy(/beam if you like).
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
#
################################################################################
import numpy as np
import sys, time, subprocess, os
sys.path.append('../utilities')
sys.path.append('../new_bethermin')
sys.path.append('../source_handling')
from genmap import genmap_gauss
import config
from clus_popmap import clus_popmap
from clus_get_clusparams import clus_get_clusparams
from clus_get_data import clus_get_data
from clus_format_bethermin import clus_format_bethermin
from astropy.io import fits
from FITS_tools.hcongrid import hcongrid
import matplotlib.pyplot as plt
import math

def clus_sim_background(genbethermin=1,fluxcut=0,saveplots=1,savemaps=0,genpowerlaw=0,\
            genradave=1,addnoise=0,yeslens=1,resolution='nr',nsim=1,bolocam=0,verbose=1,\
            errmsg=None,superplot=0):

    # clusters = ['a0370','a1689','a1835','a2218','a2219','a2390',
              # 'cl0024','ms0451','ms1054','ms1358','rxj0152','rxj1347']
    clusters = ['rxj1347']

    nclust = len(clusters)

    nsim = 100
    # fluxcut == [[10**(-np.arrage(6)+2))],0] # was used for doing radial average maps

    # Welcome the user
    print(' ')
    if verbose:
        'Welcome to background simulator v 1.0 Python Edition'
    print(' ')

    for iclust in range(nclust):

        # Set up some path handling
        bethfile = str('bethermin_'+ clusters[iclust])
        bethcat = str('bethermin_'+ clusters[iclust] + '_cat.fits')

        # First get the cluster parameters
        if verbose:
            print('Fetching cluster %s parameters' %(clusters[iclust]))

        # Use the version of clusparams that we have already worked with
        params, err = clus_get_clusparams(clusters[iclust],iclust,verbose=verbose)
        if err:
            if verbose:
                print('clus_get_clusparams exited with error: ' + err)
            exit()

        # now get the spire maps
        if verbose:
            print('Fetching SPIRE maps for %s' %(clusters[iclust]))

        maps, err = clus_get_data(clusters[iclust], iclust, resolution=resolution, bolocam=bolocam, verbose=verbose)
        if err:
            if verbose:
                print('clus_get_data exited with error: ' + err)
            exit()

        # set some initial params, these don't usually change
        ncols =  len(maps)
        wave = [250.,350.,500.] # In units of microns
        fwhm = [17.6, 23.9, 35.2]
        map_size = [300,220,150] # set by Conley sim maps generation
        bands = ncols * [0]

        if resolution == 'fr':
            '''no idea if this is still correct for new pixsize gaussian generation'''
            pixsize = np.repeat(2.,ncols)

        else:
            pixsize = [6.0, 8.33333, 12.0]

        if bolocam:
            'this is going to be different for bolocam, but not sure what...'
            pixsize = [pixsize,6]
            wave = [wave,2100.]

        for icol in range(ncols):
            if icol < 3:
                bands[icol] = maps[icol]['band']
            else:
                bands[icol] = 'BOLOCAM'

        for isim in range(nsim):
            print('Starting sim %s of %s sims' %(isim,nsim))

            if genbethermin: # use CONLEY
                if verbose:
                    print('Staring Bethermin.')

                gm = genmap_gauss(wave=wave, pixsize=pixsize, fwhm=fwhm)
                sim_maps = gm.generate(0.25,verbose=True)

            else : # use SIDES
                PLW = np.load(config.HOME + 'sides_sims/sides_PLW_sim%s.npy' %(isim),allow_pickle=True)
                PMW = np.load(config.HOME + 'sides_sims/sides_PMW_sim%s.npy' %(isim),allow_pickle=True)
                PSW = np.load(config.HOME + 'sides_sims/sides_PSW_sim%s.npy' %(isim),allow_pickle=True)
                sim_maps = [PSW,PMW,PLW]

            for icol in range(ncols):
                retcat,truthtable = clus_format_bethermin(icol,sim_maps,maps,map_size[icol],bands[icol],clusters[iclust],
                                        pixsize[0],fwhm[icol],fluxcut=fluxcut,zzero=params['z'],superplot=superplot,savemaps=savemaps,genbethermin=genbethermin)

                np.save(config.HOME + 'new_bethermin/sides_cat/sides_cat_%s_%s.npy'%(bands[icol],isim),truthtable,allow_pickle=True)

                # if yeslens == 1:
                #     print('Starting ', clusters[iclust], ' lens run for ', bands[icol], '...')
                #
                #     # create output data files for lenstool
                #     lnfile = config.HOME + 'model/' + clusters[iclust] + '/' + clusters[iclust] + '_cat.cat'
                #     ltfile = config.SIMBOX + clusters[iclust] + '_image_' + bands[icol] + '.dat'
                #
                #     # # creates a symbolic link to this file in the current working directory
                #     subprocess.Popen(['ln -s %s' %(lnfile)],shell=True)
                # 
                #     # Important Note: Lenstool will not process any sources with arcsecond positional coordinates over ~1800
                #     # # calling lenstool program, waiting to continue until program finishes
                #     subprocess.call(['/usr/local/bin/lenstool %s -n' %(config.HOME + 'model/' + clusters[iclust] + '/bestopt.par')],shell=True)
                #
                #     # this moves the output from lenstool to the directory SIMBOX and renames it
                #     subprocess.call(['mv -f image.all %s' %(ltfile)],shell=True)
                #
                #     # post process cleanup, lenstool makes these files
                #     # os.remove(config.SIM + 'mu.fits')
                #     # os.remove(config.SIM + 'image.dat')
                #     # os.remove(config.SIM + 'dist.dat')
                #     # os.remove(config.SIM + 'pot.dat')
                #     # os.remove(config.SIM + 'clump.dat')
                #     # os.remove(config.SIM + 'para.out')
                #     # os.remove(config.SIM + 'sipos.dat')
                #     # os.remove(config.SIM + 'source.dat')
                #     # os.remove(config.SIM + 'sort.dat')
                #
                # else :
                #     ltfile = config.SIMBOX + clusters[iclust] + '_image_' + bands[icol] + '.dat'
                #
                # # populate the sim maps with sources from lenstool
                # outmap,table = clus_popmap(ltfile,maps[icol],map_size[icol],bands[icol],clusters[iclust],pixsize[icol],fwhm[icol],loz=retcat,superplot=superplot,savemaps=savemaps)
                # np.save(config.HOME + 'sides_cat_%s_%s.npy'%(bands[icol],isim),table,allow_pickle=True)
                # # modify the outmap to remove pixels that have been flagged
                # # if icol < 3 :
                # #     whpl = np.where(maps[icol]['flag'] > 0) # check that this is right
                # #     outmap[whpl] = np.nan
                #
                # # reshape map to original SPIRE size
                # outmap = outmap[0:maps[icol]['signal'].shape[0],0:maps[icol]['signal'].shape[1]]
                # maps[icol]['signal'] = outmap
                #
                # # adding random noise to the signal map
                # if addnoise == 1 :
                #     print('Adding noise to sim maps')
                #     np.random.seed(102793)
                #     error =  np.array(maps[icol]['error']).flatten()
                #     signal = np.array(maps[icol]['signal']).flatten()
                #     noise = np.random.normal(loc=0.0,scale=1.0,size=(maps[icol]['signal'].shape[0],maps[icol]['signal'].shape[1])).flatten()
                #     normalization = 2*np.pi*((3)/2.355)**2
                #
                #     for i in range(len(signal)):
                #         if not math.isnan(error[i]) : #Jy/beam
                #             noise[i] = error[i]*noise[i]
                #             signal[i] = signal[i]*normalization
                #
                #         else :
                #             signal[i] = np.nan
                #             noise[i] = np.nan
                #
                #     flat = noise + signal
                #     maps[icol]['signal'] = flat.reshape(maps[icol]['signal'].shape[0],maps[icol]['signal'].shape[1])
                #
                #     # save the current sim map in /data/sim_clusters
                #     print('Savemaps set, saving maps in %s' %(config.CLUSSIMS))
                #     savefile = config.CLUSSIMS + clusters[iclust] + '/' + clusters[iclust] + '_' + bands[icol] + '_sim03%02d.fits' %(isim)
                #     # create image HDU
                #     hda = fits.PrimaryHDU(maps[icol]['signal'],maps[icol]['shead'])
                #     hdul = fits.HDUList(hdus=hda)
                #     hds = fits.ImageHDU(maps[icol]['signal'],maps[icol]['shead'])
                #     hdn = fits.ImageHDU(maps[icol]['error'],maps[icol]['shead'])
                #     hde = fits.ImageHDU(maps[icol]['exp'].data,maps[icol]['shead'])
                #     hdm = fits.ImageHDU(maps[icol]['mask'],maps[icol]['shead'])
                #     hdul.append(hds)
                #     hdul.append(hdn)
                #     hdul.append(hde)
                #     hdul.append(hdm)
                #     hdul.writeto(savefile,overwrite=True)

if __name__ == '__main__' :
    clus_sim_background(genbethermin=0,fluxcut=0,saveplots=0,savemaps=0,genpowerlaw=0,\
                        genradave=1,addnoise=1,yeslens=1,resolution='nr',nsim=0,bolocam=0,\
                        verbose=0,errmsg=None,superplot=0)
