################################################################################
# NAME : compute_rings.py
# DATE STARTED : June 24, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE : computes the radial average 
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
from math import *
import numpy as np

def compute_rings(maps, params, binwidth, superplot=0, verbose=1, noconfusion=None):
    ncols = len(maps)
    radave = np.empty(ncols)
    if noconfusion:
        confusionnoise = np.empty(ncols)
    else:
        confusionnoise = np.array([5.8, 6.3, 6.8]) * 0.001 #not sure if this will work

    bands = ['PSW', 'PMW', 'PLW']
    nbins = 0
    for i in range(ncols):
        mapsize = maps[i]['astr']['NAXIS']
        pixsize = maps[i]['pixsize']

        maxrad = ceil(pixsize * np.amax(mapsize) / sqrt(2.0))

        nbinsp = maxrad / binwidth + 1

        if nbinsp > nbins:# PURPOSE : idk lol

            nbins = nbinsp

    #!P.MULTI = [0,1,ncols] I have absolutley no idea what this means...
    # This is a function in idl that allows multiple plots on a page
    # If we want to have the plots come up on the screen as we plot them we will need
    # find a matplot replacement.
    for m in range(ncols):
        if verbose:
            print('Average for ' + maps[m]['band'])

        mapsize = maps[m]['astr']['NAXIS']
        pixsize = maps[m]['pixsize']

        maxrad = ceil(pixsize * np.amax(mapsize) / sqrt(2.0))

        radbin = []
        for i in range(int(nbins)):
            radbin.append(i / (int(nbins)-1) + 3)
        radbin[0] = 0.0
        radbin = np.array(radbin)

        #midbin = abs(TS_DIFF(radbin,1)/ 2.0) i don't understand what TS_DIFF is
        midbinp = np.empty(int(nbins))
        midwei = np.empty(int(nbins))
        for i in range(len(radbin)):
            midbin[i] = radbin[i] + midbin[i]

        fluxbin = np.empty(nbins)
        hitbin = np.empty(nbins)
        errbin = np.empty(nbins)

        #convert RA/DEC to xy coordinates

        conv = pixsize / binwidth
        calfac = (1*10**-6) / (1.13*25.0**2*( pi /180.0)**2/3600.0**2)
        confnoise = confusionnoise[m]
        mask = make_noise_mask(maps, m)
        whplerr = []
        for i in range(mask.shape[0]):
            for j in range(mask.shape[1]):
                if np.isfinite(mask[i,j]) == False:
                    whplerr.append([i,j])
                else:
                    pass
        for value in whplerr:
            mask[value] = 1
        tempmap = np.empty(mapsize[0], mapsize[1])

        for i in range(mapsize[0]):
            for j in range(mapsize[1]):
                #thisrad = REFORM(args) call to a function called REFORM, i don't know what the equivalent of REFORM would be in python because i don't understand what reform is doing in idl
                thisrad = thisrad[0]
                rad = []
                for k in range(thisrad.shape[0]):
                    for h in range(thisrad.shape[1]):
                        if abs(thisrad[k,h] - midbin[k,h]) == floor(abs(thisrad[k,h] -midbin[k,h])):
                            rad.append([k,h])
                for vlaue in rad:
                    midbinp[rad] = midbinp[rad] + thisrad #this may not work...might have to use np.add()
                    midwei[rad] = midwei[rad] + 1

                maps[m]['mask'][i,j] = maps[m]['mask'][i,j] + mask[i,j]
                tempmap[i,j] = thisrad #i don't get this.

                if len(rad) > 0 and maps[m]['mask'][i,j] == 0 and rad < maxrad:
                    thiserr = calfac * sqrt(maps[m]['error'][i,j]**2 + confnoise**2)
                    for value in rad:
                        fluxbin[value] = fluxbin[value] + calfac * maps[m]['srcrm'][i,j] / thiserr**2
                        hitbin[value] = hitbin[value] + 1.0 / thiserr**2
        clusname = maps[m]['name']
        file = '../SPIRE/cluster_analysis/plots/radmap_' + bands[m] + '_' + clusname + '.fits'
        writefits(file, data=tempmap, header_dict=maps[m]['shead'])

        whpl = []
        nanitout = []
        for i in range(midwei.shape[0]):
            for j in range(midwei.shape[1]):
                if midwei[i,j] >= 1.0:
                    whpl.append([i,j])
                else:
                    nanitout.append([i,j])

        for value in whpl:
            midbinp[value] = midbinp[value] / midwei[value]

        whpl = []
        nnan = []
        for i in range(hitbin.shape[0]):
            for j in range(hitbin.shape[1]):
                if hitbin >= 0:
                    whpl.append([i,j])
                else:
                    nnan.append([i,j])

        for value in whpl:
            fluxbin[value] = fluxbin[value] / hitbin[value]
            errbin[value] = sqrt(1.0 / hitbin[value])\
        if len(nnan) > 0:
            for value in whnan:
                midbinp[value] = np.nan
                fluxbin[value] = np.nan
                errbin[value] = np.nan
        radave[m] = {'band' : maps[m]['band'],
                     'midbin' : midbinp,
                     'fluxbin' : fluxbin,
                     'errbin' : errbin}
        if superplot:
            #function to plot the things can probably use matplotlib here.

    return radave
