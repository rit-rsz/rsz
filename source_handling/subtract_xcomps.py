################################################################################
# NAME : subtract_xcomps.py
# DATE STARTED : June 21, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE : idk lol
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
from math import *

def subtract_xcomps(maps, simflag=0, verbose=1):
    ncols = len(maps)

    if maps[0]['band'] != 'PSW':
        err = 'first element of map structure is not PSW, aborting.'
        if verbose:
            print(err)
        return None, err

    for i in range(ncols):
        if verbose:
            print('On band %s' %(maps[i]['band']))
        width = sqrt(maps[icol]['widtha']**2 - maps[0]['widtha']**2) / maps[0]['pixsize']
        #kern = PSF_GAUSSIAN we need to have a function for PSF_GAUSSIAN
        inmap = maps[0]['srcrm']
        maps[0]['xclean'] = maps[0]['srcrm']
        whpl = []
        whnan = []
        for j in range(inmap.shape[0]):
            for k in range(inmap.shape[1]):
                if np.isfinite(inmap[j,k]) == True:
                    whpl.append([j,k])
                else:
                    whnan.append([j,k])
        for indexes in whnan:
            maps[0]['mask'][indexes] = 1
            inmap[indexes] = 0.0
        #xmap = CONVOL which is some math thing that i don't understand what it does so i can't find an equivalent.
        #then it does a where in xmap and calls HASTROM.
        #meanclip is called which supposedly just returns a mean, but i don't know if its got some weird thing in it.
        whpl = []
        whnan = []
        for j in range(maps[i]['srcrm'].shape[0]):
            for k in range(maps[i]['srcrm'].shape[1]):
                if np.isfinite(maps[icol]['srcrm'][j,k]) == True and np.isfinite(xmap[j,k]) == True:
                    whpl.append([j,k])
                else:
                    whnan.append([j,k])
        #then it calls plot, so i don't know if we can just matplotlib here but ok i guess

        #coeff = Least squares fit of some other stuff.

        x = np.empty(1000)
        for j in range(x.shape[0]):
            x[j] = (j - 500.0) / 1000.0
        #another call to plot this new thing stil ldon't know if we can just use matplotlib here.
        for j in range(maps[i]['xclean'].shape[0]):
            for k in range(maps[i]['xclean'].shape[1]):
                maps[i]['xclean'][j,k] = maps[i]['srcrm'][j,k] - (coeff[1] * xmap_align[j,k] + coeff[0])

        for indexes in whnan:
            maps[i]['mask'][indexes] = 1

        datasub = maps[i]['xclean']

        if not simflag:
            filename = config.CLUSDATA + 'sz/' + maps[i]['name'] + str(maps[i]['band']) + '_xc.fits'
        else:
            filename = config.CLUSDATA + 'sz/sim/' + maps[i]['name'] + str(maps[i]['band']) + '_xc.fits'
        #we do writefits here again so i really need to find an equivalent...

        #another call to contour need to fix this.

    for i in range(maps[0]['xclean'].shape[0]):
        for j in range(maps[0]['xclean'].shape[1]):
            maps[0]['xclean'][i,j] = maps[0]['xclean'][i,j] - mean(maps[0]['xclean']) # this is questionable for a 2d array.

    return maps, err
