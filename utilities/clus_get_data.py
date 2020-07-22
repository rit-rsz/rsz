
################################################################################
# NAME : get_data.py
# DATE STARTED : June 11, 2019
# AUTHORS : Dale Mercado & Benjamin Vaughan
# PURPOSE : This is where the map data is retrieved from various cluster catalogs
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
#   Victoria Butler - 7/19 - edits to calfac for unit conversion
################################################################################
import numpy as np
import matplotlib.pyplot as plt
from math import *
from astropy.io import fits
import os
import sys
from get_spire_beam import *
from get_spire_beam_fwhm import *
import config
import matplotlib.pyplot as plt
sys.path.append('/home/bjv7945/rsz/Planck_Model')
from cirrus_priors import one_map

def clus_get_data(clusname, isim,manpath=0, resolution = 'nr', bolocam=None,
            verbose = 1, version = '1', manidentifier=None, sgen=None):

    errmsg = False
    #if the manual path flag is not set then go the normal route.
    if manpath == 0 and sgen is None:
        hermfiles = []
        hermdir = config.CLUSDATA + 'hermes_clusters/'
        hermcount = 0
        #check in the HerMES directory.
        for x in os.listdir(hermdir):
            if x.startswith(clusname):
                if resolution in x and version in x:
                    hermfiles.append(hermdir + x) #this SHOULD work
                    hermcount += 1
        #check in the HLS directory.
        hlsdir = config.CLUSDATA + 'hls_clusters/'
        hlsfiles = []
        hlscount = 0
        for x in os.listdir(hlsdir):
            if x.startswith(clusname):
                if x.startswith(clusname):
                    hlsfiles.append(hlsdir + x)
                    hlscount += 1
        #check in the snapshot directory.
        snapdir = config.CLUSDATA + 'snapshot_clusters/'
        snapfiles = []
        snapcount = 0
        for x in os.listdir(snapdir):
            if x.startswith(clusname):
                snapfiles.append(snapdir + x)
                snapcount += 1

        #if it did not find a PSw, PMW, and PLW file for the given cluster in any of these directories then return an error.
        if snapcount ^ hermcount ^ hlscount != 3 :
            errmsg = ('Problem finding exclusive files, file dump is: %s %s %s' %
                        (snapfiles, hermfiles, hlsfiles))
            if verbose:
                print(errmsg)
            return None, errmsg


        elif hermcount == 3:
            files = hermfiles
            nfiles = hermcount
            #we are using HerMES data.
            print('hermes data')

        elif hlscount == 3:
            files = hlsfiles
            nfiles = hlscount
            #we are using HLS data.
            print('hls data')

        elif snapcount == 3:
            files = snapfiles
            nfiles = snapcount
            #we are using snapshot data.
            print('snap data')

    #this is in the case there is a manual path that needs to be specified. (currently I don't think we need this functionality.)
    elif manpath == 1:
        mancount = 0
        manfiles = []
        file_flag = False
        for x in os.listdir(manpath):
            if clusname in x and manidentifier in x:
                manfiles.append(manpath + x)

        if len(manfiles) != 3:
            errmsg = 'Cannot find files in %s' %(manpath)
            if verbose:
                print(errmsg)
            return None, errmsg
        nfiles = len(manfiles)
        files = manfiles

    elif sgen is not None:
        #this is for looking for simulation data.
        '''
        sgen =  1: old bethermin (bad)
                2: old bethermin (good)
                3: SIDES or Conley
        '''
        simfiles = []
        simdir = config.CLUSSIMS + clusname + '/'
        simcount = 0
        for x in os.listdir(simdir):
            if x.startswith(clusname):
                if ('0' + str(sgen)+ str(isim).zfill(2)) in x and 'fits' in x and 'BOLOCAM' not in x:
                    simfiles.append(simdir + x)
                    simcount += 1
        files = simfiles
        nfiles = simcount

    if bolocam:
        nfiles = nfiles +1

    maps = []

    for i in range(nfiles):
        if i < 3:
            maps.append(clus_read_file(files[i], clusname, verbose=verbose,sgen=sgen))

        else:
            maps[ifile] = np.empty(clus_read_bolocam(clusname,verbose=verbose)) #args need to be filled in bolocam one


    #the purpose of the below code is to organize the maps objects so that our program doesn't bork out and think the PSW
    #map is the PLW map or the PLW map is the PMW map, etc.
    sort_order = {'PSW' : 0, 'PMW' : 1, 'PLW' : 2}

    maps.sort(key = lambda x : sort_order[x['band']])

    return maps, errmsg

##################################################################################################
##################################################################################################

def clus_read_file(file, clusname, verbose=0, sgen=None):


    #set the band based off of what is in the filename
    if 'PSW' in file:
        band = 'PSW'
    elif 'PMW' in file:
        band = 'PMW'
    elif 'PLW' in file:
        band = 'PLW'

    #collect data from the files file.
    hdul = fits.open(file)

    if sgen == None :
        img = hdul[1] #image object
        err = hdul[2] #error map
        exp = hdul[3] #exposure map
        mask = hdul[4] #mask map
        #only do the below steps for real data.
        planck, iris = one_map(name, map)
        planck_hdu = fits.ImageHdu(planck, shead)
        iris_hdu   = fits.ImageHdu(iris, shead)
        hdul.append(planck_hdu)
        hdul.append(iris_hdu)
        hdu.writeto(file, overwrite=True)
    elif int(sgen) == 3 :
        img = hdul['signal'] #image object
        err = hdul['error'] #error map
        exp = hdul['exp'] #exposure map
        mask = hdul['mask'] #mask map
        noise = hdul['noise']

    else :
        img = hdul[1] #image object
        err = hdul[2] #error map
        exp = hdul[3] #exposure map
        mask = hdul[4] #mask map

    #this adds noise into our map if it's a simulation. (I think we are currently doing this elsewhere.)
    if sgen == 1 or sgen == 2: # if using the old pipeline bethermin sims
        mapsize = img.data.shape
        noisemap = 1.0 * err.data * np.random.standard_normal((mapsize[0], mapsize[1]))
        img.data = img.data + noisemap
        img.data = img.data - np.nanmean(img.data)
    if sgen == 3:
        img.data = img.data - np.nanmean(img.data)

    if 'CDELT1' in img.header.keys():
        #calculating the pixsize based off of astrometry parameters.
        pixsize = 3600 * mean([abs(img.header['CDELT1']),abs(img.header['CDELT2'])])
        #if the given astrometry has CDELT values and not a cd matrix create a cd matrix.
        img.header['cd_1'] = img.header['CDELT1']
        img.header['cd_2'] = 0
        img.header['cd_1'] = 0
        img.header['cd_2'] = img.header['CDELT2']
        #don't think we need to do this for the error header as well.
        # err.header['cd_1'] = err.header['CDELT1']
        # err.header['cd_2'] = 0
        # err.header['cd_1'] = 0
        # err.header['cd_2'] = err.header['CDELT2']

    else:
        #if the astrometry is given with a cd matrix calculate the pixsize this way.
        pixsize = 3600 * \
                  mean([abs(img.header['CD1_1']+img.header['CD2_1']), \
                        abs(img.header['CD2_1'] + img.header['CD2_2'])])

    psf = get_spire_beam(pixsize=pixsize, band=band, verbose=0) #creates a 2D gaussian as our psf.
    widtha = get_spire_beam_fwhm(band) #arcsecs (fwhm of gaussian)
    width = widtha / (sqrt(8 * log(2)) / pixsize) # width in pixels #sqrt(8 * log(2)) is conversion to sigma of gauss
    calfac = 1 / (config.calfac * (get_spire_beam_fwhm(band))**2) #calibration factor based off of FWHM of our beam.
#   Gets header information from a fits image. Astropy should be able to do this
    #This is potentially depreciated.


    head = img.header
    herr = err.header

    srcrm = np.zeros(img.shape)
    xclean = np.zeros(img.shape)
    flag = np.zeros(img.shape)

    #put everything into a dictionary.
    maps = {'name':clusname,
          'file':file,
          'band':band,
          'signal':img.data,
          'srcrm':srcrm,
          'xclean':xclean,
          'iris': iris,
          'planck': planck,
          'error':err.data,
          'mask':mask.data,
          'noise' : noise.data,
          'flag':flag,
          'exp':exp,
          'shead':head,
          'ehead':herr,
          'pixsize':pixsize,
          'psf':psf,
          'width':width,
          'widtha':widtha,
          'calfac':calfac,
          'JY2MJy':config.JY2MJy}


    return maps


# This can be done with numpy
def mean(data):
    length = len(data)
    sum = 0
    for i in range(length):
        sum += data[i]
    average = sum / length
    return average

if __name__ == '__main__':
    clus_get_data('a0370', nsim='197')
