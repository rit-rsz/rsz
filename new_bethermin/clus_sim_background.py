
################################################################################
# NAME : clus_sim_background.py
# DATE STARTED : September 27, 2019
# AUTHORS : Dale Mercado
# PURPOSE : The wrapper for Conley's Bethermin simmulation generator. fluxcut is in Jy(/beam if you like).
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
import sys
sys.path.append('../utilities')
import config
from clus_popmap import clus_popmap
from clus_get_clusparams import clus_get_clusparams
from clus_get_data import clus_get_data
from clus_format_bethermin import clus_format_bethermin

def clus_sim_background(clusname,genbethermin=1,fluxcut=0,saveplots=1,savemaps=0,genpowerlaw=0,\
            genradave=1,addnoise=0,yeslens=1,resolution='fr',nsim=0,bolocam=0,verbose=1,errmsg=None):

    clusters = ['a0370','a1689','a1835','a2218','a2219','a2390',
              'cl0024','ms0451','ms1054','ms1358','rxj0152','rxj1347']

    clusters = clusname #Is this how you select the one that you want to be working on???

    nclust = len(clusters)

    nsim = 300
    fluxcut == [[10**(-np.arrage(6)+2))],0]

    for iclust in range(nclust):

        # Line for the output save file path

        for isim in range(nsim):

            print('On sim ' + str(isim+1) + 'of' + str(nsims))

            if icut == 0:
                genbethermin == 1
            else:
                genbethermin == 0

        # Welcome the user
        print(' ')
        if verbose:
            'Welcome to background simulator v 1.0 Python Edition'
        print(' ')

        # Set up ps plotting if appropriate
        if saveplots:
            '''Need to create Nuplot for this aswell'''

        # Set up some path handling
        bethfile = str('bethermin_'+ clusname)
        bethcat = str('bethermin_'+ clusname + '_cat.fits')

        # First get the cluster parameters
        if verbose:
            print('Fetching cluster parameters')
        # Use the version of clusparams that we have already worked with
        params, err = clus_get_clusparams(clusname,verbose=verbose)
        if err:
            if verbose:
                print('clus_get_clusparams exited with error: ' + err)
            exit()

        # now get the spire maps
        if verbose:
            print('Fetching SPIRE maps.')
        maps, err = clus_get_data(clusname=clusname, resolution=resolution, bolocam=bolocam, verbose=verbose)
        if err:
            if verbose:
                print('clus_get_data exited with error: ' + err)
            exit()


        ncols =  len(maps)
        wave = [250.,350.,500.] # In units of um
        fwhm = [17.6, 23.9, 35.2]

        if resolution == 'fr':
            '''no idea if this is still correct for new pixsize gaussian generation'''
            pixsize = np.repeat(2.,ncols)

        else:
            pixsize = [6.0, 8.33333, 12.0]

        if bolocam:
            'this is going to be different for bolocam, but not sure what...'
            pixsize = [pixsize,6]
            wave = [wave,2100.]

        for icol in range(ncol):
            if icol < 3:
                bands[icol] = maps['band']
            else:
                bands[icol] = 'BOLOCAM'

            psffiles[icol] = str(config.CLUSBOX + bands[icol] + '_psf.fits')

        if genbethermin:
            if verbose:
                print('Setting up bethermin model.')

            print('Staring Bethermin.')

            sim_maps = genmap(wave,pixsize,fwhm)

            ncols == 1
            bethcat = catfile
            bethfile = outfile


            for icol in range(ncols):
                clus_format_bethermin(icol,STRING(!CLUSSBOX,bethcat),STRING(!CLUSSBOX,bethfile),
                    bands[icol],clusname,pixsize[icol],FLUXCUT=fluxcut,ZZERO=params.z,RETCAT=lozcat)

                if yeslens == 1:

                    print('Starting ', clusname, ' lens run for ', bands[icol], '...')

                    # create output data files for lenstool
                    lnfile = config.CLUSHOME + 'model/' + clusname + '/' + clusname + '_cat.cat'
                    ltfile = config.CLUSSBOX + clusname + '_image_' + bands[icol] + '.dat'

                    # creates a symbolic link to this file in the current working directory
                    subprocess.Popen(['ln -s %s' %(lnfile)],shell=True)

                    # calling lenstool program
                    subprocess.Popen(['/usr/local/bin/lenstool', config.CLUSHOME + 'cluster_analysis/model/' + clusname + '/bestopt.par', '-n'],shell=True)

                    # i'm assuming that this moves the output from lenstool to the directory CLUSSBOX and renames it
                    subprocess.Popen(['mv image.all %s' %(ltfile)],shell=True)

                    # post process cleanup
                    os.remove(clusname+'_cat.cat')
                    os.remove('*.fits')
                    os.remove('*.dat')
                    os.remove('*.out')

                else :
                    ltfile = config.CLUSHOME + 'model/' + clusname + '/' + clusname + '_cat.cat'

                # populate the sim maps with sources from lenstool
                clus_popmap(ltfile,maps[icol],MAP=outmap,LOZ=lozcat)

                if icol < 3 :
                    whpl = np.where(maps[icol]['flag'] > 0) # check that this is right
                    outmap[whpl] = np.nan

                maps[icol]['signal'] = outmap


                # adding random noise to the signal map
                if addnoise == 1 :
                    print('Adding noise to sim maps')

                    confnoise = 0.0
                    noisemap = np.random.random(seed,) #noisemap = RANDOMN(SEED,(*maps[icol]).astr.naxis)
                    (maps[icol]['signal'] = maps[icol]['signal'] + (maps[icol]['mask']*maps[icol]['error'])*noisemap)

                if savemaps == 1 :
                    # save the current sim map in /data/sim_clusters

                elif savemaps == 2 :
                    # save the current sim map in the 0200 naming convention in /data/bethermin_sims

                if saveplots == 1 :
                    # make some plots
