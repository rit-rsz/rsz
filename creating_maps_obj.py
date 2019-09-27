#HOW TO CALL:
#import this script.
# maps = create_map_obj(filename, dir)
#filename is a list of all three files
# i.e. [PSW_file, PMW_file, PLW_file]
#directory is the directory where all these files are kept
#example is at the bottom.



import scipy.io
import numpy as np
# from config import * #(this line will give me access to all directory variables)
import matplotlib.pyplot as plt
from math import *
from astropy.io import fits
import scipy.signal
import os
# from astropy.modeling import models, fitting
from astropy.convolution import Gaussian2DKernel, Gaussian1DKernel


def rebin(a, new_shape):
    shape = a.shape
    M = int(shape[0])
    N = int(shape[1])
    m, n = new_shape
    print('M', M, 'N', N, 'm', m, 'n', n)
    print('M/m', M/m, 'N/n', N/n)
    if m<M:
        return a.reshape((m,int(M/m),n,int(N/n))).mean(3).mean(1)
    else:
        return np.repeat(np.repeat(a, m/M, axis=0), n/N, axis=1)


def get_spire_beam_fwhm(band):
    if band == 'PSW':
        beamFWHM = 18.0
        return beamFWHM
    elif band == 'PMW':
        beamFWHM = 25.0
        return beamFWHM
    elif band == 'PLW':
        beamFWHM = 36.0
        return beamFWHM
    elif band == 'BOLOCAM':
        beamFWHM = 144.0
        return beamFWHM
    else:
        print('Unknown band:', band)

def get_spire_beam(band=None, pixsize=0,npixx=0, npixy=0,
                   xcent=0, ycent=0,bolometer=0, fwhm='',
                   norm=0, oversamp=0, verbose=1,
                   factor=0):
    errmsg = False
    if bolometer != 0:
        if verbose:
            print('PSFs for specific bolometer not yet supported -- ignoring')

    # Check if we have been given a band
    # If not make an assumption
    if band == None:
        if verbose:
            print('Band parameter not supplied, assuming PSW')
        band = 'PSW'
        # Check if we have been givin a pixel size
    if pixsize == 0:
        # band = upper(band)
        # units arcsec/pixel
        if band == 'PSW':
            pixsize = 6
        if band == 'PMW':
            pixsize = 8 + (1/3)
        if band == 'PLW':
            pixsize = 12
        else:
            print('Unknown band'+band)
        if verbose:
            print('pixsize paramter not supplied assuming %s arcsec' , (pixsize))

    # Figure out which color this is
    # Units should be in arcsec
    if len(fwhm) == 0:
        beamFWHM = get_spire_beam_fwhm(band)
    else:
        beamFWHM = fwhm

    # This seems to be redundent but is in the idl script
    if beamFWHM < 0:
        print('Invalid Beam FWHM value'+ str(beamFWHM))

    # Check if we've been given the map size, if not assume something
    # npixx/npixy will be the final number of pixels
    if npixx == 0:
        npixx = round(beamFWHM * 5 / (pixsize))
        if npixx % 2 != 1:
            npixx +=1
            if verbose:
                print('npixx not supplied, using ",I0"')
    #If no y size then assume same as x
    if npixy == 0:
        npixy = npixx

    #Make sure that these have been cast from a float properly or we get errors
    npixx = int(ceil(npixx))
    npixy = int(ceil(npixy))

    if npixx % 2 != 1 and verbose:
        print('WARNING: npixx not odd, so PSF will not be centered')
    if npixy % 2 != 1 and verbose:
        print('WARNING: npixy not odd, so psf will not be centered')

    # Now deal with oversampling
    if oversamp == 0:
        ioversamp = 7
    else:
        ioversamp = round(oversamp) #in case user provides float
        if ioversamp % 2 != 1:
            print('Oversamp must be an odd interger!')

    x_gen = npixx * ioversamp
    y_gen = npixy * ioversamp
    gen_pixsize = np.float64(pixsize) / ioversamp

    # Check if we have been givin the center, if not assume middle
    if xcent == 0:
        if verbose:
            print('xcent parameter not supplied, assuming array center')
        ixcent = x_gen / 2
    else:
        # Adjust for oversampling
        ixcent = xcent * ioversamp
        if ioversamp > 1:
            ixcent = ixcent + ioversamp / 2

    if ycent == 0:
        if verbose:
            print('ycent parameter not supplied, assuming array center')
        iycent = y_gen / 2
    else:
        iycent = ycent * ioversamp
        if ioversamp > 1:
            iycent = iycent +ioversamp / 2

    # Normalize FWHM to pixels
    beamFWHM /= gen_pixsize
    # Convert the FWHM to a standard deviation for astropy fit. From psf_gaussian
    stdev = beamFWHM / (sqrt(8 * log(2)))

    # If we want this normalized then call with norm flag set
    if factor:
        # 1D beam
        beamkernraw = Gaussian1DKernel(stdev, x_size=x_gen)
        beamkern = np.array(beamkernraw)
        if norm:
            beamkern = beamkern / beamkern.max()
        if ioversamp > 1:
            beamkern = rebin(beamkern,(npixx,npixy))
    else:
        beamkernraw = Gaussian2DKernel(stdev,x_size = x_gen, y_size =  y_gen)
        beamkern = np.array(beamkernraw)
        if norm:
            beamkern = beamkern / beamkern.max()
        if ioversamp > 1:
            beamkern = 	rebin(beamkern,(npixx,npixy))



    # print('RIGHT HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!', beamkern)
    # # Use for debugging
    # # plt.plot for 1d
    # # plt.plot(beamkern, drawstyle='steps')
    # # plt.imshow for 2d & colorbar
    # plt.imshow(beamkern, interpolation='none', origin='lower')
    # plt.colorbar()
    # plt.xlabel('x [pixels]')
    # plt.ylabel('y [pixels]')
    # plt.show()
    beamkern = 1

    return beamkern

def mean(data):
    length = len(data)
    sum = 0
    for i in range(length):
        sum += data[i]
    average = sum / length
    return average

def clus_read_file(file,band,clusname,verbose=0,simmap=0):
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
    map = hdul[1]
    err = hdul[2]
    exp = hdul[3]
    flag = hdul[4]

    if simmap:
        mapsize = map.data.shape
        noisemap = 1.0 * err.data * np.random.standard_normal((mapsize[0], mapsize[1]))
        map.data = map.data + noisemap
        map.data = map.data - np.nanmean(map.data)
        # maps.append(thismap['thismap'])

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
    calfac = 1 / (50* (get_spire_beam_fwhm(band))**2)
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
          'error':err.data, #check
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
          'JY2MJy':50} #check


    return maps

def create_map_obj(filname, dir):
    band = ['PSW', 'PMW', 'PLW']
    maps = [[],[],[]]
    for i in range(3):
        maps[i] = clus_read_file(dir + filename[i], band[i], 'some cluster name')
#band name should be 'PSW, PMW, PLW'
    return maps
#clusname doesn't matter unless your script cares about clusname in which cause put
#the name of the cluster you are using,


if __name__ == '__main__':
    dir = '/data/mercado/SPIRE/bethermin_sims/a0370/' #change this to whatever your directory for keeping files is
    filename = ['a0370_PSW_sim0110.fits', 'a0370_PMW_sim0110.fits', 'a0370_PLW_sim0110.fits'] #this is a list of the filenames.
#please put the file names in order of PSW, PMW, PLW
    print(create_map_obj(filename, dir))
