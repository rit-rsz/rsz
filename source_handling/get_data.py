
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
import scipy.io
import numpy as np
# from config import * #(this line will give me access to all directory variables)
import matplotlib.pyplot as plt
from math import *
from astropy.io import fits
import os
import sys
#import pyfits
sys.path.append('../utilities')
from get_spire_beam import *
from get_spire_beam_fwhm import *
import config
import matplotlib.pyplot as plt

def get_data(clusname, manpath=0, resolution = 'nr', bolocam=None,
            verbose = 1, version = '1', manidentifier=None, nsim=0):
    # place holders for the if statements to work until I add them to the input for get data
    # This will happen when the script is fully functional
    print('clus_get_data')
    errmsg = False
    if not bolocam:
        cols = ['PSW','PMW','PLW']
        bolocam = 0
    else:
        cols = ['PSW','PMW','PLW','BOLOCAM']
        bolocam = 1

#   If there is no manual path set
    if manpath == 0 and sim == 0:
        hermfiles = []
        hermdir = config.CLUSDATA + 'hermes_clusters/'
        hermcount = 0
        for x in os.listdir(hermdir):
            if x.startswith(clusname):
                if resolution in x and version in x:
                    hermfiles.append(hermdir + x) #this SHOULD work
                    hermcount += 1

        hlsdir = config.CLUSDATA + 'hls_clusters/'
        hlsfiles = []
        hlscount = 0
        for x in os.listdir(hlsdir):
            if x.startswith(clusname):
                if x.startswith(clusname):
                    hlsfiles.append(hlsdir + x)
                    hlscount += 1

        snapdir = config.CLUSDATA + 'snapshot_clusters/'
        snapfiles = []
        snapcount = 0
        for x in os.listdir(snapdir):
            if x.startswith(clusname):
                snapfiles.append(snapdir + x)
                snapcount += 1

        if snapcount ^ hermcount ^ hlscount != 3 :
            errmsg = ('Problem finding exclusive files, file dump is:', \
                        snapfiles, hermfiles, hlsfiles)
            if verbose:
                print(errmsg)
            return None, errmsg
        if hermcount == 3:
            files = hermfiles
            nfiles = hermcount
            print('hermes data')

        if hlscount == 3:
            files = hlsfiles
            nfiles = hlscount
            print('hls data')

        if snapcount == 3:
            files = snapfiles
            nfiles = snapcount
            print('snap data')

    #   Manual Path option
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


    if nsim != 0:
        simfiles = []
        simdir = config.CLUSSIMS + clusname + '/'
        simcount = 0
        for x in os.listdir(simdir):
            if x.startswith(clusname):
                if nsim in x and 'fits' in x and 'BOLOCAM' not in x:
                    simfiles.append(simdir + x)
                    simcount += 1
        files = simfiles
        nfiles = simcount

    if bolocam:
        nfiles = nfiles +1

    print(simfiles)
    maps = []


#   Need to tweek the syntax of this for loop
    for ifile in range(nfiles):
        if ifile < 3: # I feel like this is supposed to be ifile <= 3 :
            if any(col in files[ifile] for col in cols): # checks if name of band is in the filename
                maps.append(read_file(files[ifile], cols[ifile], clusname, verbose=verbose))

            else:
                errmsg = 'Problem finding ' + cols[ifile] + ' file.'
                if verbose:
                    print(errmsg)
        else:
                maps[ifile] = np.empty(clus_read_bolocam(clusname,verbose=verbose)) #args need to be filled in bolocam one

    new_file = 'None'

    print(len(maps))
    for i in range(len(maps)):
        print(maps[i]['file'])
        if 'PSW' in maps[i]['file'] and i != 0:
            new_file = maps[0]['file']
            holder = maps[0]
            maps[0] = maps[i]
            maps[0]['band'] = 'PSW'


    if 'PMW' in new_file:
        holder2 = maps[1]
        maps[1] = holder
        maps[1]['band'] = 'PMW'
        new_file = holder2['file']

    elif 'PLW' in new_file:
        holder2 = maps[2]
        maps[2] = holder
        maps[2]['band'] = 'PLW'
        new_file = holder2['file']

    if 'PLW' in new_file:
        maps[2] = holder2
        maps[2]['band'] = 'PLW'


    elif 'PMW' in new_file:
        maps[1] = holder2
        maps[1]['band'] = 'PMW'

    for i in range(len(maps)):
        if 'PMW' in maps[i]['file']:
            new_file = maps[1]['file']
            holder = maps[1]
            maps[1] = maps[i]
            maps[1]['band'] = 'PMW'

    if 'PLW' in new_file:
        maps[2] = holder
        maps[2]['band'] = 'PLW'


    return maps, errmsg

##################################################################################################
##################################################################################################

def read_file(file,band,clusname,verbose=0):
    # Should this calfac be listed as a self.calfac or not?
    '''
    Calfac has been added to config.py as a constant.
    This is the first place it is created a used.
    '''
    # calfac = (pi/180.0) * (1/3600.0)**2 * (pi / (4.0 * log(2.0))) * (1e6)
    # This will ultimatly be in the list of constants
    # The rest of the scrpit involves idl_libs stuff that
    # will get grabbed from astropy
    hdul = fits.open(file)
    # print(test.header)
    map = hdul[1]
    # print(map.header)
    # exit()
    err = hdul[2]
    exp = hdul[3]
    flag = hdul[4]

    if 'CDELT1' in map.header.keys():
        pixsize = 3600 * mean([abs(map.header['CDELT1']),abs(map.header['CDELT2'])])
        map.header['cd_1'] = map.header['CDELT1']
        map.header['cd_2'] = 0
        map.header['cd_1'] = 0
        map.header['cd_2'] = map.header['CDELT2']
        # err.header['cd_1'] = err.header['CDELT1']
        # err.header['cd_2'] = 0
        # err.header['cd_1'] = 0
        # err.header['cd_2'] = err.header['CDELT2']

    else:
        pixsize = 3600 * \
                  mean([abs(map.header['CD1_1']+map.header['CD2_1']), \
                        abs(map.header['CD2_1'] + map.header['CD2_2'])])

    psf = get_spire_beam(pixsize=pixsize, band=band)
    #psf = 4 #for xid test only...
    widtha = get_spire_beam_fwhm(band) #arcsecs (sigma of gaussian)
    width = widtha / (sqrt(8 * log(2)) * pixsize) # width in pixels
    calfac = 1 / (config.calfac * (get_spire_beam_fwhm(band))**2)
#   This should be defined in the catsrsc file
    '''
    This was moved to config.py
    '''
    # JY2MJy = 1e6


#   Gets header information from a fits image. Astropy should be able to do this
    astr = {}
    try:
        cd11 = map.header['CD1_1']
        cd12 = map.header['CD1_2']
        cd21 = map.header['CD2_1']
        cd22 = map.header['CD2_2']
    except KeyError:
        pass # i don't like the way this is coded probably have to change it later
    for keys in map.header.keys():
        if 'NAXIS' in keys:
            astr.update({keys : map.shape})
        if 'CD1_1' in keys:
            x =  np.array([[cd11, cd12], [cd21, cd22]])
            astr.update({'CD' : x})
        if 'CDELT' in keys:
            astr.update({keys : map.header[keys]})
        if 'CRPIX1' in keys:
            x = np.array([map.header['CRPIX1'], map.header['CRPIX2']])
            astr.update({'CRPIX' : x})
        if 'CTYPE1' in keys:
            x = np.array([map.header['CTYPE1'], map.header['CTYPE2']])
            astr.update({'CTYPE' : x})
        if 'CRVAL1' in keys:
            x = np.array([map.header['CRVAL1'], map.header['CRVAL2']])
            astr.update({'CRVAL' : x})
        if 'LONGPOLE' in keys:
            astr.update({keys : map.header[keys]})
        if 'LATPOLE' in keys:
            astr.update({keys : map.header[keys]})


    head = map.header
    herr = err.header

    print(type(head))
#   Not sure if this is the correct syntax for astr naxis
    srcrm = np.zeros(astr['NAXIS'])
    xclean = np.zeros(astr['NAXIS'])
    mask = np.zeros(astr['NAXIS'])

#  Need to generate default mask map
#  whnan = WHERE(FINITE(map) EQ 0,countnan)
#  IF countnan GT 0 THEN mask[whnan] = 1

    maps = {'name':clusname, #check
          'file':file, #check
          'band':band, #check
          'signal':map.data, #check
          'srcrm':srcrm, #check
          'xclean':xclean, #check
          'error':err, #check
          'mask':mask, #check
          'flag':flag, #check
          'exp':exp, #check
          'shead':head, #nope
          'ehead':herr, #nope
          'astr':astr, #check
          'pixsize':pixsize, #check
          'psf':psf, #check    dumb = fits.PrimaryHDU(maps['signal'], map.header)
          'width':width, #check
          'widtha':widtha, #check
          'calfac':calfac, #check
          'JY2MJy':config.JY2MJy} #check

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
    get_data('a0370', nsim='197')
