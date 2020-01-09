
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
import math
import sys, time, subprocess, os
sys.path.append('../utilities')
sys.path.append('../new_bethermin')
sys.path.append('../source_handling')
sys.path.append('../sz')
from genmap import genmap_gauss
import config
from clus_popmap import clus_popmap
from clus_get_clusparams import clus_get_clusparams
from clus_get_data import clus_get_data
from clus_format_bethermin import clus_format_bethermin
from astropy.io import fits
from FITS_tools.hcongrid import hcongrid
from astropy.convolution import convolve, Gaussian2DKernel, Tophat2DKernel
from astropy.modeling.models import Gaussian2D
from clus_add_sziso import clus_add_sziso
import matplotlib.pyplot as plt

def clus_sim_background(genbethermin=1,fluxcut=0,saveplots=1,savemaps=0,genpowerlaw=0,\
            genradave=1,addnoise=0,yeslens=1,resolution='nr',nsim=1,bolocam=0,verbose=1,\
            errmsg=None,superplot=0):

    clusters = ['a0370']#,'a1689','a1835','a2218','a2219','a2390',
             # 'cl0024','ms0451','ms1054','ms1358','rxj0152','rxj1347']

    nclust = len(clusters)
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

        print('made it out of clus_get_data')
        for isim in range(nsim):
            print('Starting sim %s of %s sims' %(isim,nsim))

            if genbethermin:
                if verbose:
                    print('Staring Bethermin.')

                gm = genmap_gauss(wave=wave, pixsize=pixsize, fwhm=fwhm)
                sim_maps = gm.generate(0.25,verbose=True)

            for icol in range(ncols):
                lozcat = clus_format_bethermin(icol,sim_maps,maps,bands[icol],clusters[iclust],
                                        pixsize[icol],fwhm[icol],fluxcut=fluxcut,zzero=params['z'],superplot=superplot,savemaps=savemaps)

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
                    # os.remove(lnfile)
                    # # creates a symbolic link to this file in the current working directory
                    subprocess.Popen(['ln -s %s' %(lnfile)],shell=True)
                    #
                    # # calling lenstool program, waiting to continue until program finishes
                    subprocess.call(['/usr/local/bin/lenstool %s -n' %(config.HOME + 'model/' + clusters[iclust] + '/bestopt.par')],shell=True)
                    #
                    # this moves the output from lenstool to the directory SIMBOX and renames it
                    subprocess.call(['mv image.all %s' %(ltfile)],shell=True)

                    # post process cleanup
                    # os.remove(config.SIM + 'mu.fits')
                    # os.remove(config.SIM + 'image.dat')
                    # os.remove(config.SIM + 'pot.dat')
                    # os.remove(config.SIM + 'dist.dat')
                    # os.remove(config.SIM + 'para.out')

                else :
                    ltfile = config.SIMBOX + clusters[iclust] + '_image_' + bands[icol] + '.dat'
                    # ltfile = config.SIMBOX + 'a0370_image_PSW_test.dat'

                print('fwhm_sim: ',pixsize[icol],fwhm[icol])
                # populate the sim maps with sources from lenstool
                outmap = clus_popmap(ltfile,maps[icol],bands[icol],clusters[iclust],pixsize[icol],fwhm[icol],loz=lozcat,superplot=superplot,savemaps=savemaps)
                # os.remove(ltfile)
                # modify the outmap to remove pixels that have been flagged
                # if icol < 3 :
                #     whpl = np.where(maps[icol]['flag'] > 0) # check that this is right
                #     outmap[whpl] = np.nan

                # reshape map to original SPIRE size
                # outmap = outmap[0:maps[icol]['signal'].shape[0],0:maps[icol]['signal'].shape[1]]
                # kern = Gaussian2DKernel(3.0)
                # outmap_smooth = convolve(outmap,kern)
                # maps[icol]['signal'] = outmap_smooth
                maps[icol]['signal'] = outmap
                print('map size:',maps[icol]['signal'].shape[0],maps[icol]['signal'].shape[1])
                plt.imshow(outmap)
                plt.savefig('sim_map_%s.png' %(bands[icol]))

                # adding random noise to the signal map
                if addnoise == 1 :
                    print('Adding noise to sim maps')
                    ''' under the assumption that the error map is exclusively instrument noise'''
                    # confnoise = [4.8,4.4,4.8] #mJy/beam  #[5.8,6.3,6.8] * 1e-3 / 2. #(SQRT(8.))
                    np.random.seed(102793)
                    error =  np.array(maps[icol]['error']).flatten()
                    signal = np.array(maps[icol]['signal']).flatten()
                    noise = np.random.normal(scale=1.0,size=(maps[icol]['signal'].shape[0],maps[icol]['signal'].shape[1])).flatten()

                    for i in range(len(signal)):
                        if not math.isnan(error[i]) : #Jy/beam
                            noise[i] = error[i]*noise[i]
                            # signal[i] = signal[i] * maps[icol]['calfac'] / 1e6
                            # signal[i] = signal[i]*(1.13*(fwhm[icol]/pixsize[icol])**2)
                            signal[i] = signal[i]*6.0
                            # if bands[icol] == 'PMW':
                            #     signal[i] = signal[i]*1.9
                            # if bands[icol] == 'PLW':
                            #     signal[i] = signal[i]*4.0
                        else :
                            signal[i] = np.nan
                            noise[i] = np.nan

                    print('lenght',len(noise),len(signal))
                    flat = noise + signal
                    noise_map = noise.reshape(maps[icol]['signal'].shape[0],maps[icol]['signal'].shape[1])
                    maps[icol]['signal'] = flat.reshape(maps[icol]['signal'].shape[0],maps[icol]['signal'].shape[1])
                    signal_map = signal.reshape(maps[icol]['signal'].shape[0],maps[icol]['signal'].shape[1])

                    if savemaps == 1 or savemaps == 2 :
                        # save the current sim map in /data/sim_clusters
                        print('Savemaps set, saving maps in %s' %(config.SIMBOX+'sim_clusters/'))
                        savefile = config.SIMBOX + 'sim_clusters/' + clusters[iclust] + '_' + bands[icol] + '_noise.fits'
                        # create image HDU
                        hda = fits.PrimaryHDU(signal_map,maps[icol]['shead'])
                        hdul = fits.HDUList(hdus=hda)
                        hds = fits.ImageHDU(maps[icol]['signal'],maps[icol]['shead'],name='s+n')
                        hdul.append(hds)
                        hdn = fits.ImageHDU(noise_map,maps[icol]['shead'],name='noise')
                        hdul.append(hdn)
                        hdul.writeto(savefile,overwrite=True)

                        # TESTING ######################################################
                        hdun = fits.PrimaryHDU(noise_map,maps[icol]['shead'])
                        hdun.writeto(config.SIMBOX+'noise_%s.fits' %(bands[icol]),overwrite=True)
                        hdun = fits.PrimaryHDU(signal_map,maps[icol]['shead'])
                        hdun.writeto(config.SIMBOX+'sim_%s.fits' %(bands[icol]),overwrite=True)
                        hdun = fits.PrimaryHDU(maps[icol]['signal'],maps[icol]['shead'])
                        hdun.writeto(config.SIMBOX+'both_%s.fits' %(bands[icol]),overwrite=True)
                        ###############################################################

                    # if savemaps == 1 :
                    #     # save the current sim map in the 0200 naming convention in /data/bethermin_sims
                    #     print('Savemaps set , saving maps in %s' %(config.CLUSNSIMS))
                    #     savefile = config.CLUSNSIMS + clusters[iclust] + '_' + bands[icol] + '_sim0' + str(isim) + '.fits'
                    #     # hdx = fits.PrimaryHDU(maps[icol]['signal'],maps[icol]['shead'])
                    #     hdx = fits.PrimaryHDU(third_map,maps[icol]['shead'])
                    #     if os.path.isfile(savefile):
                    #         os.remove(savefile)
                    #     hdx.writeto(savefile,overwrite=True)

                    # if saveplots == 1 :
                    #     make some plots
                    # maps[icol]['signal'] = fits.getdata(config.SIMBOX+'both_%s.fits' %(bands[icol]))
                    # signal = maps[icol]['signal'].flatten()
                    # for i in range(len(signal)):
                    #     if bands[icol] == 'PMW':
                    #         signal[i] = signal[i]*1.9
                    #     if bands[icol] == 'PLW':
                    #         signal[i] = signal[i]*4.0
                    # maps[icol]['signal'] = signal.reshape(maps[icol]['signal'].shape[0],maps[icol]['signal'].shape[1])

            # yin_coeff = [2.50,1.91,2.26,3.99,1.36,2.42,1.59,1.90,3.99]
            # yin = [x*1e-4 for x in yin_coeff]
            # tin = [7.2,10.1,7.7,9.8,4.5,8.6,7.8,5.5,10.9]
            # new_maps,junk = clus_add_sziso(maps,yin=yin,tin=tin)
            # print(len(new_maps))
            # hdun = fits.PrimaryHDU(new_maps[0]['signal'],maps[0]['shead'])
            # hdun.writeto(config.SIMBOX+'sz_%s.fits' %(bands[0]),overwrite=True)
            # hdun = fits.PrimaryHDU(new_maps[1]['signal'],maps[1]['shead'])
            # hdun.writeto(config.SIMBOX+'sz_%s.fits' %(bands[1]),overwrite=True)
            # hdun = fits.PrimaryHDU(new_maps[2]['signal'],maps[2]['shead'])
            # hdun.writeto(config.SIMBOX+'sz_%s.fits' %(bands[2]),overwrite=True)
            exit()

if __name__ == '__main__' :
    clus_sim_background(genbethermin=1,fluxcut=0,saveplots=1,savemaps=1,genpowerlaw=0,\
                        genradave=1,addnoise=1,yeslens=0,resolution='nr',nsim=25,bolocam=0,\
                        verbose=0,errmsg=None,superplot=0)
