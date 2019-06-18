import numpy as np
from astropy.io import fits
import random
import scipy.io as sp

#fits.open() may need to be replaced or expanded on depending on how things go


def get_simmaps(simflag=1, verbose=1, sb=0, xc=0, clusname, nsim): # clusname, SIMFLAG=simflag. NSIM=nsim, SB=sb, XC=xc,
#success = success, errmsg = errmsg, verbose = verbose
    #verbose = 1 in args
    #simflag = 1 in args
    #sb = 0 in args
    #xc = 0 in args
    if simflag > 0 and len(nsim): #what is n_elements?
        print('simmap set but nsim not supplied! Aborting.')
    if simflag == 0:
        nsim = np.nan
    cols = ['PSW', 'PMW', 'PLW']
    bigpix = [6.0,8.3,12.0]
    ncols = len(cols)
    maps = np.empty(ncols)
    if simflag == 1:
        for i in range(ncols-1):
            mapfilename = clusdata + 'sim_clusters/' + clusname + '_' + cols[i] +
            '.sav'
            thismap = sp.readsav(mapfilename)
            print(thismap.calcfac) #where does thismap come from?
            thismap.calcfac = 1.0 / (thismap.calcfac * (thismap.pixsize)*
            thismap.pixsize)
            if sb:
                overmap = fits.open(clusdata + 'sim_clusters/' + clusname +
                '_' + cols[i] + 'sb.fits', head) #read fits needs to be replaced
                thismap.xclean = overmap
                whpl = np.where(np.isfinite(overmap) == False)
                for index in whpl:
                    thismap.mask[index] = 1
            elif xc:
                overmap = fits.open(clusdata + 'sim_clusters/' + clusname + '_' +
                cols[i] + 'xc.fits', head) #read fits needs to be replaced
                thismap.xclean = overmap
                whpl = np.where(np.isfinite(overmap) == False)
                for index in whpl:
                    thismap.mask[index] = 1
            maps[i] = thismap #maybe needs PTR_NEW equivalent maybe not?
    else:
        for i in range(n_cols-1):
            thissim = 'nsim' #idk what format = I04 does...
            mapfilename = clusdata + 'bethermin_sims/' + clusname +'/' + clusname +
            '_' + cols[i] + '_sim' + thissim + '.sav'
            thismap = sp.readsav(mapfilename)
            mapsize = thismap.signal.shape #still don't know what thismap is.
            noisemap = 1.0 * thismap.error * np.random.standard_normal((mapsize[0],
            mapsize[1]))
            thismap.signal = thismap.signal + noisemap #STILL don't know what thismap is
            maps[i] = thismap
    return maps
