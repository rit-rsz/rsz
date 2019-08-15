################################################################################
# NAME : compute_rings.py
# DATE STARTED : June 24, 2019
# AUTHORS : Benjamin Vaughan & Victoria Butler
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
from astropy import wcs
from astropy.io import fits

def compute_rings(maps, params, binwidth, superplot=0, verbose=1, noconfusion=None):

    # init params
    bands = ['PSW', 'PMW', 'PLW']
    ncols = len(maps)
    radave = np.empty(ncols)
    if noconfusion:
        confusionnoise = np.empty(ncols)
    else:
        ar = [5.8, 6.3, 6.8]
        confusionnoise = np.array([x*0.001 for x in ar])

    # create bins for different metrics (mid,flux,err,rad,hit,max)
    for m in range(ncols):
        if verbose:
            print('Average for ' + maps[m]['band'])

        # make sizes for different radial averages
        mapsize = maps[m]['astr']['NAXIS']
        pixsize = maps[m]['pixsize']
        maxrad = ceil(pixsize * np.amax(mapsize) / sqrt(2.0))
        nbinsp = maxrad / binwidth + 1

        # make array objects to fill later
        midbinp = np.empty(len(nbinsp))
        midwei = np.empty(len(nbinsp))
        fluxbin = np.empty(len(nbinsp))
        hitbin = np.empty(len(nbinsp))
        errbin = np.empty(len(nbinsp))

        # make radbin
        stop = maxrad * (len(nbinsp)/(nbinsp-1)) + 3
        radbin = np.linspace(0.0,stop,len(nbinsp))

        # make midbin
        midbin = np.absolute([x / 2.0 for x in np.diff(radbin,1)])
        midbin = np.add(radbin,midbin)

        ''' This whole section is gonna need some serious TLC'''
        #convert RA/DEC to pixel coordinates
        ''' Need to check that fidrad and fidded are in degrees'''
        ra = params['fidrad'] * u.deg
        dec = params['fidded'] * u.deg
        c = SkyCoord(ra, dec)
        # grabbing header from maps file
        hdul = fits.open(maps[m]['file'])
        w = wcs(hdul[1].header)
        #converting ra/dec to pixel coords
        px, py = skycoord_to_pixel(c, w)
        # ===============================================================


        conv = pixsize / binwidth
        ''' Not sure if this is the right calfac'''
        calfac = (1*10**-6) / (1.13*25.0**2*( pi /180.0)**2/3600.0**2)
        confnoise = confusionnoise[m]
        ''' We haven't made this script yet...
            Nick made a version but not sure if it works'''
        mask = make_noise_mask(maps, m)
        # ===============================================================

        # use mask to identify bad pixels
        for i in range(mask.shape[0]):
            for j in range(mask.shape[1]):
                if np.isfinite(mask[i,j]) == False:
                    mask[i,j] = 1
        # ====================================================

        # not exactly sure what this is doing but looks like radius of bin rings
        tempmap = np.empty(mapsize[0], mapsize[1])
        for i in range(mapsize[0]):
            for j in range(mapsize[1]):
                nthisrad = 0 # keep a counter for how many times we find a minimum
                maps[m]['mask'][i,j] = maps[m]['mask'][i,j] + mask[i,j]
                thisrad = pixsize * np.sqrt(np.add((np.power((i - px),2),np.power((j - py),2))))
                #thisrad = REFORM(args) # removes the first dimension of the array (1,10,10) --> (10,10)
                thisrad = thisrad[0] # why do we care about only the first value ?
                tempmap[i,j] = thisrad
                for k in range(len(midbin)):
                    if abs(thisrad[k] - midbin[k]) == min(abs(thisrad[k] -midbin[k])):
                        nthisrad += 1
                        midbinp[k] = midbinp[k] + thisrad
                        midwei[k] = midwei[k] + 1
                        if nthisrad != 0 and maps[m]['mask'][i,j] == 0 and (k <= maxrad) :
                            thiserr = calfac * sqrt(maps[m]['error'][i,j]**2 + confnoise**2)
                            fluxbin[k] = fluxbin[k] + calfac * maps[m]['srcrm'][i,j] / thiserr**2
                            hitbin[k] = hitbin[k] + 1.0 / thiserr**2

# ------------------------------------------------------------------------------------------------------------

        clusname = maps[m]['name']
        file = '../SPIRE/cluster_analysis/plots/radmap_' + bands[m] + '_' + clusname + '.fits'
        writefits(file, data=tempmap, header_dict=maps[m]['shead'])

        for i in range(len(nbinsp)):
            if midwei[i] >= 1.0 :
                midbinp[i] = midbinp[i] / midwei[i]
            if hitbin[i] >= 0 :
                fluxbin[i] = fluxbin[i] / hitbin[i]
                errbin[i] = sqrt(1.0 / hitbin[i])

        if len(nnan) > 0:
            for value in whnan:
                midbinp[value] = np.nan
                fluxbin[value] = np.nan
                errbin[value] = np.nan



        # save new bin data to dictionary & return to fitsz
        radave[m] = {'band' : maps[m]['band'],
                     'midbin' : midbinp,
                     'fluxbin' : fluxbin,
                     'errbin' : errbin}
        if superplot:
            plt.plot(radave[m]['midbin'],)

    return radave

if __name__ == '__main__' :
    maps,err = get_data('a2218')
    # params gets used as input for AD2XY later
    compute_rings(maps,params)
