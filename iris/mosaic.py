################################################################################
# Name : mosaiq
# Purpose : copy of the IDL routine mosaique_iris.pro ported over to python
# this code does a transformation on an image to rotate and align it with
# another image.
#Author : Benjamin Vaughan
#Start Date : September 23, 2019
#Additional Info
#
################################################################################
from astropy.io import fits
import numpy as np
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
from astropy.wcs import WCS as world
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK4
from astropy.coordinates import Galactic
from astropy.coordinates import BarycentricTrueEcliptic
import astropy.units as u
import sys
from math import *
import os
from bilinear import mbilinear
from utility import *
from get_iris_map import get_iris
import matplotlib.pyplot as plt
import iris_config

def mosaic(header, band=4, catname=iris_config.IrisLookupFile, dir=iris_config.IrisDir):
    '''
    purpose: create mosaic of images
    Inputs : Header - header information -should contain some astrometry information
             band   - the band you wish to look at
                - 1 : 12 micron
                - 2 : 25 micron
                - 3 : 60 micron
                - 4 : 100 micron (default)
             catname- filename of the catalog
             dir    - directory where iris maps are
    Outputs:
    '''


    header.set('NAXIS', 2)
    w = world(header)
    try:
        equinox = header['EQUINOX']
    except KeyError:
        equinox = header['EPOCH']


    x_size = header['NAXIS1']
    y_size = header['NAXIS2']
    print('Handling %s elements' % (x_size * y_size))
    x = np.arange(0, x_size)
    y = np.arange(0, y_size)
    xlist = np.tile(1, x_size)
    ylist = np.tile(1, y_size)
    xmap = np.outer(x, ylist)
    ymap = np.outer(xlist, y)

    result = np.zeros((x_size, y_size))
    weight = np.zeros((x_size, y_size))
    new_c = pixel_to_skycoord(xmap, ymap, w, origin=0)
    #this is converting the pixel coords to right ascension and declination in fk4
    ra = np.asarray(new_c.ra.to_string(decimal=True), dtype=float)
    dec = np.asarray(new_c.dec.to_string(decimal=True), dtype=float)

    ctype = get_cord_type(header)


    ra = nan2undef(ra)
    dec = nan2undef(dec)
    ra = np.asarray(ra)
    dec = np.asarray(dec)
    #converted these to arrays for easier data manipulation

    inum, ramin, ramax, raavg, decmin, decmax, decavg, medianval, noise_key = np.loadtxt(catname, unpack=True)
    numel = int(inum[-1]) #number of entries in the text file
    print('Checking for ISSA maps that intersect with the given header')

    ind1 = np.where(ra != -32768)[0]
    ind2 = np.where(dec != -32768)[0]

    ind = []
    for i in range(ra.shape[0]):
        for j in range(ra.shape[1]):
            if ra[i,j] != -32768 and dec[i,j] != -32768:
                ind.append([i,j])
    ind = np.asarray(ind)
    good_inds = np.zeros(numel)

    c1min = np.min(ra[ind])
    c1max = np.max(ra[ind])
    c2min = np.min(dec[ind])
    c2max = np.max(dec[ind])

    for i in range(numel):
        if c1min > ramin[i] and c1min < ramax[i] and c2min > decmin[i] and c2min < decmax[i]:
            good_inds[i] = 1
        elif c1min > ramin[i] and c1min < ramax[i] and c2max > decmin[i] and c2max < decmax[i]:
            good_inds[i] = 1
        elif c1max > ramin[i] and c1max < ramax[i] and c2max > decmin[i] and c2max < decmax[i]:
            good_inds[i] = 1
        elif c1max > ramin[i] and c1max < ramax[i] and c2min > decmin[i] and c2min < decmax[i]:
            good_inds[i] = 1
        elif ramin[i] > c1min and ramin[i] < c1max and decmin[i] > c2min and decmin[i] < c2min:
            good_inds[i] = 1
        elif ramax[i] > c1min and ramax[i] < c1max and decmin[i] > c2min and decmin[i] < c2min:
            good_inds[i] = 1
        elif ramin[i] > c1min and ramin[i] < c1max and decmax[i] > c2min and decmax[i] < c2min:
            good_inds[i] = 1
        elif ramax[i] > c1min and ramax[i] < c1max and decmax[i] > c2min and decmax[i] < c2min:
            good_inds[i] = 1

    good_inds = np.where(good_inds > 0)[0]
    if good_inds.shape[0] == 0:
        print('No ISSA map corresponds to the header given')
        exit()

    print('%s ISSA maps will be combined to produce the mosaic' %(good_inds.shape[0]))

    for i in range(good_inds.shape[0]):
        mapi = get_iris(inum[good_inds[i]], dir=dir, band=band)


        try:
            del(mapi[0].header['NAXIS3'])
        except KeyError:
            pass

        #do the transform back to pixel coords
        w = world(mapi[0].header)
        x, y = skycoord_to_pixel(new_c, w, origin=0)
        tempo = mbilinear(x, y, mapi[0].data)
        for j in range(tempo.shape[0]):
            for k in range(tempo.shape[1]):
                if tempo[j,k] != -32768:
                    weight[j,k] = weight[j,k] + 1
                    result[j,k] = result[j,k] + tempo[j,k]
    indw = []
    complement = []
    for i in range(weight.shape[0]):
        for j in range(weight.shape[1]):
            if weight[i,j] > 0:
                result[i,j] = result[i,j] / weight[i,j]
            else:
                result[i,j] = -32768
    #in interpolating the image the x and y axes get flipped.
    result = np.swapaxes(result, 0, 1)
    if band == 4:
        result -= 0.65 #subtract CIB to get just dgl / cirrus


    return result
