################################################################################
# NAME : subtract_cat.py
# DATE STARTED : June 21, 2019
# AUTHORS : Benjamin Vaughan & Victoria Butler
# PURPOSE : the purpose of this script is to subtract the data from XID from
#           the original signal
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
#   VLB - 8/13/19 - trimmed down the code to get rid of unneccessary for loops,
#                   added contour plots to visualize output, robust deletion and
#                   creation of new fits files for debugging
#   VLB - 10/15/19 - changing to work with PCAT_SPIRE residual maps instead of XID+
################################################################################
import numpy as np
import sys
import matplotlib.pyplot as plt
# from clus_get_data import *
from astropy.io import fits
sys.path.append('/utilities')
import config
sys.path.append(config.HOME + 'multiband_pcat')
from pcat_spire import lion
from save_fits import writefits
import os, time
import multiprocessing as mp
from multiprocessing import Pool
from clus_make_noise_mask import clus_make_noise_mask

def clus_subtract_cat(maps, dI, verbose=1, nsim=0, saveplot=0, superplot=0):
    err = None

    ''' Multithreading PCAT Setup
    # default is to have all of the 3 maps returned
    # manager = mp.Manager()
    # ret_maps = manager.dict()
    # resid_maps = []
    # a = time.time()
    # p1 = mp.Process(target=run_pcat, args=(maps,ret_maps))
    # resid_maps.append(p1)
    # p2 = mp.Process(target=run_pcat, args=(maps,ret_maps))
    # p3 = mp.Process(target=run_pcat, args=(maps,ret_maps))
    # p4 = mp.Process(target=run_pcat, args=(maps,ret_maps))
    # p1.start()
    # p2.start()
    # p3.start()
    # p4.start()
    # p1.join()
    # p2.join()
    # p3.join()
    # p4.join()
    # b = time.time()
    # print('###########################')
    # print('TIME ELAPSED:' , b-a)
    # print('###########################')
    '''

    ### This is for saving residuals from pcat so we can do testing without having to rerun through pcat.
    # plt.imsave('/home/butler/rsz/pcat_resid_0.png',resid_maps[0],format='png')
    # plt.imsave('/home/butler/rsz/pcat_resid_1.png',resid_maps[1],format='png')
    # plt.imsave('/home/butler/rsz/pcat_resid_2.png',resid_maps[2],format='png')
    # np.save('/home/butler/rsz/pcat_resid.npy',resid_maps)

    ###This is for loading in the saved maps so we can run catsrc without running pcat.
    # map1 = np.load('/home/butler/rsz/pcat_resid.npy',allow_pickle=True)
    # map2 = np.load('/home/butler/rsz/pcat_resid1.npy',allow_pickle=True)
    # map3 = np.load('/home/butler/rsz/pcat_resid2.npy',allow_pickle=True)
    # resid_maps = [map1[0],map2[0],map3[0]]

    #make a noise mask and populate error map with nans where mask is 1 so they don't effect the fit.
    for i in range(len(maps)):
        mask = clus_make_noise_mask(maps, i)
        maps[i]['mask'] = maps[i]['mask'] + mask
        for j in range(maps[i]['signal'].shape[0]):
            for k in range(maps[i]['signal'].shape[1]):
                if maps[i]['mask'][j,k] == 1:
                    maps[i]['error'][j,k] = 0.0

    plt.imshow(maps[i]['error'])
    plt.colorbar()
    plt.savefig('checking_pcat.png')
    plt.clf()

    resid_maps = run_pcat(maps)

    # make sure both input maps exist
    # if maps.any() == None or resid_maps.any() == None :
    #     return None, 'clus or pcat_spire map structure input is absent'

    for i in range(len(resid_maps)): # only loop through how many residuals we have
        if verbose==1:
            print('Setting Subtracted Maps for %s' %(maps[i]['band']))

        # make the difference image
        datasub = resid_maps[i]
        # make full data map object
        datafull = maps[i]['signal']


        if superplot:
            # plotting for debugging purposes
            fig1, ax1 = plt.subplots()
            fig1, ax2 = plt.subplots()
            cp1 = ax1.contour(datafull)
            ax1.set_title('clus_subtract_cat: Signal map for %s' %(maps[i]['band']))
            cp2 = ax2.contour(datasub)
            ax2.set_title('clus_subtract_cat: Catalog subtracted map for %s' %(maps[i]['band']))
            plt.show()

            plt.imshow(datasub)
            plt.colorbar()
            plt.title('clus_subtract_cat: Catalog subtracted map for %s' %(maps[i]['band']))
            plt.show()

            # save new subtracted signal to map
            # reshape pcat square map to original SPIRE size
        maps[i]['srcrm'] = datasub[0:maps[i]['signal'].shape[0],0:maps[i]['signal'].shape[1]]

        ''' Testing '''
        filename = config.HOME + 'outputs/pcat_residuals/' + maps[i]['name'] + '_' + maps[i]['band'] + '_residx25' + '.fits'
        hda = fits.PrimaryHDU(maps[i]['srcrm'],maps[i]['shead'])
        hda.writeto(filename,overwrite=True)
        plt.imshow(maps[i]['srcrm'])
        plt.savefig('after_pcat_%s.png'%(maps[i]['band']))
        plt.clf()

        # if saveplot:
        #     if nsim != 0:
        #         filename = config.HOME + 'outputs/pcat_residuals/' + maps[i]['name'] + '_' + maps[i]['band'] + '_' + str(nsim) + '.fits'
        #     else:
        #         filename = config.HOME + 'outputs/pcat_residuals/' + maps[i]['name'] + '_' + maps[i]['band'] + '_' + '.fits'
        #
        #     hda = fits.PrimaryHDU(maps[i]['srcrm'],maps[i]['shead'])
        #     hda.writeto(filename,overwrite=True)


    return maps, err

# def run_pcat(maps,ret_maps):
def run_pcat(maps):
    ob = lion(band0=0, band1=1, band2=2, map_object=maps, auto_resize=True, make_post_plots=True, openblas=True, cblas=False, nsamp=10, residual_samples=1, visual=False)
    resid_maps = ob.main()
    # ret_maps[0] = resid_maps
    return resid_maps

if __name__ == '__main__' :
    import sys
    sys.path.append('../utilities')
    sys.path.append('../source_handling')
    from clus_get_data import clus_get_data
    maps, err = clus_get_data(clusname='a0370')
    resid_maps = fits.open('../new_bethermin/lens_model/nonlensedmap_a0370_PSW2.fits')[0].data
    resid_maps = resid_maps[0:263, 0:243]
    maps[0]['signal'] = resid_maps
    clus_subtract_cat(maps)
