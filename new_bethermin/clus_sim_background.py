
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

def clus_sim_background(genbethermin=1,fluxcut=0,saveplots=1,savemaps=0,genpowerlaw=0,\
            genradave=1,addnoise=0,yeslens=1,resolution='fr',nsim=1,bolocam=0,verbose=1,errmsg=None):

    clusters = ['a0370','a1689','a1835','a2218','a2219','a2390',
              'cl0024','ms0451','ms1054','ms1358','rxj0152','rxj1347']

    nclust = len(clusters)

    nsim = 300
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
        params, err = clus_get_clusparams(clusters[iclust],verbose=verbose)
        if err:
            if verbose:
                print('clus_get_clusparams exited with error: ' + err)
            exit()

        # now get the spire maps
        if verbose:
            print('Fetching SPIRE maps for %s' %(clusters[iclust]))

        maps, err = clus_get_data(clusname=clusters[iclust], resolution=resolution, bolocam=bolocam, verbose=verbose)
        if err:
            if verbose:
                print('clus_get_data exited with error: ' + err)
            exit()

        # set some initial params, these don't usually change
        ncols =  len(maps)
        wave = [250.,350.,500.] # In units of microns
        fwhm = [17.6, 23.9, 35.2]
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

            if genbethermin:
                if verbose:
                    print('Staring Bethermin.')

                gm = genmap_gauss(wave=wave, pixsize=pixsize, fwhm=fwhm)
                sim_maps = gm.generate(0.25,verbose=True)

            for icol in range(ncols):
                lozcat = clus_format_bethermin(icol,sim_maps,maps,bands[icol],clusters[iclust],
                                        pixsize[icol],fluxcut=fluxcut,zzero=params['z'])

                if yeslens == 1:
                    ''' format bethermin needs to spit out only one cat file for the band it's on
                        then you need to loop through lenstool and make it do a new run for each catfile.
                        lenstool is dumb and doesn't know which band it's on, it just re-reads the same
                        file each time, so we need to make sure we re-write which one we are on
                    '''
                    print('Starting ', clusters[iclust], ' lens run for ', bands[icol], '...')

                    # create output data files for lenstool
                    lnfile = config.HOME + 'model/' + clusters[iclust] + '/' + clusters[iclust] + '_cat.cat'
                    ltfile = config.SIMBOX + clusters[iclust] + '_image_' + bands[icol] + '.dat'
                    #
                    #
                    # # creates a symbolic link to this file in the current working directory
                    subprocess.Popen(['ln -s %s' %(lnfile)],shell=True)
                    #
                    # # calling lenstool program, waiting to continue until program finishes
                    subprocess.call(['/usr/local/bin/lenstool %s -n' %(config.HOME + 'model/' + clusters[iclust] + '/bestopt.par')],shell=True)
                    #
                    # # i'm assuming that this moves the output from lenstool to the directory CLUSSBOX and renames it
                    subprocess.call(['mv image.all %s' %(ltfile)],shell=True)

                    # post process cleanup
                    ''' deletes lenstool output files from current working directory (/home/person/rsz/new_bethermin)'''
                    # os.remove(config.SIM + clusters[iclust] + '_cat.cat')
                    # os.remove(config.SIM + '*.fits')
                    # os.remove(config.SIM + '*.dat')
                    # os.remove(config.SIM + '*.out')

                else :
                    ltfile = config.SIM + 'model/' + clusters[iclust] + '/' + clusters[iclust] + '_cat.cat'

                # populate the sim maps with sources from lenstool
                outmap = clus_popmap(ltfile,maps[icol],bands[icol],clusters[iclust],pixsize[icol],loz=lozcat)
                #
                # if icol < 3 :
                #     whpl = np.where(maps[icol]['flag'] > 0) # check that this is right
                #     outmap[whpl] = np.nan
                #
                maps[icol]['signal'] = outmap


                # adding random noise to the signal map
                # if addnoise == 1 :
                #     print('Adding noise to sim maps')
                #
                #     confnoise = 0.0
                #     noisemap = np.random.random(seed,) #noisemap = RANDOMN(SEED,(*maps[icol]).astr.naxis)
                #     (maps[icol]['signal'] = maps[icol]['signal'] + (maps[icol]['mask']*maps[icol]['error'])*noisemap)

                # if savemaps == 1 :
                #     # save the current sim map in /data/sim_clusters
                #
                # elif savemaps == 2 :
                #     # save the current sim map in the 0200 naming convention in /data/bethermin_sims
                #
                # if saveplots == 1 :
                #     # make some plots

if __name__ == '__main__' :
    clus_sim_background(genbethermin=1,fluxcut=0,saveplots=1,savemaps=0,genpowerlaw=0,\
                        genradave=1,addnoise=0,yeslens=1,resolution='nr',nsim=1,bolocam=0,\
                        verbose=1,errmsg=None)
