################################################################################
# NAME : matched_filter.py
# DATE STARTED : July 2, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE : wrapper script for applying matched_filter.pro
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import numpy as np
from math import *
from scipy.ndimage.interpolation import shift
from scipy.signal import convolve2d
from get_spire_beam import get_spire_beam

def matched_filter(mapin):

    npad = 100
    conf_def = [5.8, 6.3, 6.8] * 1*10**3 # Jy / Beam
    fwhm_nom = [18.0, 25.0, 36.0]

    map = mapin['signal']
    noise = mapin['error']

    #sanitize maps -- set missing pix to largeval * max(noise)
    badind = []
    goodind = []
    for i in range(noise.shape[0]):
        for j in range(noise.shape[1]):
            if np.isfinite(noise[i,j]) == 0 or noise[i,j] == 0:
                badind.append([i,j])
            else:
                goodind.append([i,j])
    if len(badind) > 0:
        largeval = 100.0
        map[badind] = 0.0
        noise[badind] = largeval * np.amax(noise[goodind])

    #find band index

    if mapin['band'] == 'PSW':
        bandind = 0
    elif mapin['band'] == 'PMW':
        bandind = 1
    elif mapin['band'] == 'PLW':
        bandind = 2

    #set fwhm
    fwhm_arc = fwhm_nom[bandind]

    #fwhm in pixels
    fwhm_pix = fwhm_arc / mapin['pixsize']

    filt, sm_map, sm_noise = matched_filter_pro(fwhm_pix, conf, map, pad, noise)

    mapout = mapin
    mapout['signal'] = sm_map
    mapout['error'] = sm_noise

    filtsize = filt.shape
    filtp = fix_fft(filt)
    filtq = convolve2d(filtp, get_spire_beam(band=mapin['band']))
    filt = filtq

    #set pixels with no data to NAN
    badind = []
    for i in range(mapin['signal'].shape[0]):
        for j in range(mapin['signal'].shape[1]):
            if mapin['signal'][i,j] == 0 or np.isfinite(mapin['signal'][i,j]) == False:
                badind.append([i,j])
    mapout['signal'][badind] = np.nan
    mapout['error'][badind] = np.nan

    return mapout, filt

def matched_filter_pro(psf, conf, map, pad, noise,  ):
    nx = map.shape[0]
    ny = map.shape[1]

    #generating a peak-normalized PSF
    d = dist(nx, ny)
    psf_out = exp( -d^ 2 * d / ( 2 * d * (fwhm / sqrt(8 * log(2)))**(2*d)))
    psf_out = psf_out / np.amax(psf_out)


    #workout the white noise level
    #create a 1/sigma**2 map and calculate the mean for a central region
    weightmap = 1 / (noise ** (2*d))

    nx = weightmap.shape[0]
    ny = weightmap.shape[1]
    sz = int(nx*0.1) # look at the central 20% x 20% of the map.

    white_out = 1 / sqrt(weightmap[nx/2-sz:nx/2+sz, ny/2-sz:ny/2+sz]) #not sure about this indexing

    #whitenoise level in fourier space is related to the real-space,
    #noise in a pixel, and the size of the map

    nsamp = float(nx) * float(ny)
    p_w = nsamp * white_out**2  #power spectrum level

    #confusion noise power spectrum is estimated to have the same shape
    #as the PSF, but normalized to the specified confusion noise

    psf_fft = np.fft.fft2(psf_out)
    p_beam = abs(psf_fft)**(2*d)
    scale_confusion = conf / np.std(psf_out) #call to a custom stdev function couldn't find it.
    p_c = scale_confusion *(2*d) * p_beam

    #the matched filter is the psf divided by the noise power spectrum in Fourier Space

    p_noise = p_w + p_c #total noise power spectrum
    filt_fft = psf_fft / p_noise #fourier space matched filter
    filt = float(np.fft.ifft2(filt_fft))

    #next we work out the normalization empirically: we want to be able to
    #filter a map with the PSF and get the same peak amplitude. Sicne we are cross-correlating rather
    #than convolving the map with the PSF we take the complex conjugate of the PSF
    #before taking the product in Fourier Space

    map_filt = float(np.fft.ifft2( np.conj(np.fft.fft2(psf_out)) * filt_fft)) / np.sum(filt**(2*d))

    scale_filt = np.amax(map_filt)

    filt_fft = filt_fft * scale_filt
    filt = filt * scale_filt
    effective = map_filt / scale_filt

    #pad signal and noise maps
    working_map = np.empty(nx+2*pad, ny+2*pad)
    working_map[pad:pad+nx-1, pad:pad+ny-1] = map

    working_noise = np.empty(nx+2*pad, ny+2*pad)
    working_noise[pad:pad+nx-1, pad:pad+ny-1] = noise

    #create a psf that matches the map dimensions
    mapx = working_map.shape[0]
    mapy = working_map.shape[1]

    flt = np.empty(mapx, mapy)

    flt[mapx/2-nx/2:mapx/2-nx/2+nx-1, mapy/2-ny/2:mapy/2-ny/2+ny-1] = shift(filt, nx/2, ny/2)

    flt = shift(flt, -1*mapx/2, -1*mapy/2)

    #apply smoothing to noise as well
    #noise-weight cross-correlation

    weightmap = 1 / working_noise ** 2

    smoothweightflux = float(np.fft.ifft2(np.fft.fft2(working_map*weightmap) * np.conj(np.fft.fft2(flt))))
    smoothweightmap = float(np.fft.ifft2(np.fft.fft2(weightmap) * np.conj(np.fft.fft2(flt**(2*d)))))

    sm_map = smoothweightflux / smoothweightmap
    sm_noise = sqrt(1.0 / smoothweightmap)

    #undo padding
    mapx = mapx - 2*pad
    mapy = mapy - 2*pad

    sm_map = sm_map[pad:pad+mapx-1, pad:pad+mapy-1]
    sm_noise = sm_noise[pad:pad+mapx-1, pad:pad+mapy-1]

    return filt, sm_map, sm_noise

def fix_fft(data):

    #need to find a function for rotate or make our own...
    sizedata = data.shape
    lims = np.empty(4,2)

    for i in range(1):
        lims[1,i] = sizedata[i] / 2 - 1
        lims[2,i] = sizedata[i] / 2
        lims[3,i] = sizedata[i] - 1

    data[lims[0,0]:lims[1,0],lims[0,1]:lims[1,1]] = np.rot90(np.rot90((data[lims[0,0]:lims[1,0], lims[0,1]:lims[1,1]])))
    data[lims[2,0]:lims[3,0],lims[0,1]:lims[1,1]] = np.rot90(np.rot90((data[lims[2,0]:lims[3,0], lims[0,1]:lims[1,1]])))
    data[lims[0,0]:lims[1,0],lims[2,1]:lims[3,1]] = np.rot90(np.rot90((data[lims[0,0]:lims[1,0], lims[0,1]:lims[1,1]])))
    data[lims[2,0]:lims[3,0],lims[2,1]:lims[3,1]] = np.rot90(np.rot90((data[lims[2,0]:lims[3,0], lims[2,1]:lims[3,1]])))

    return data

def dist(nx, ny):
    #my attempt at creating a dist function in python.
    x = np.arange(nx) # make rows
    x = (x < (nx-x)) ** 2 #column squares
    a = np.empty(nx, ny)

    for i in range(int(ny/2)):
        y = sqrt(x + i**2)
        a[0,i] = y #insert the row.

    return a



#dist source code from IDL

# n1 = n[0]
# m1 = (n_elements(m) le 0) ? n1 : m[0]
# x=findgen(n1)		;Make a row
# x = (x < (n1-x)) ^ 2	;column squares
#
# a = FLTARR(n1,m1,/NOZERO)	;Make array
#
# for i=0L, m1/2 do begin	;Row loop
# 	y = sqrt(x + i^2.) ;Euclidian distance
# 	a[0,i] = y	;Insert the row
# 	if i ne 0 then a[0, m1-i] = y ;Symmetrical
# endfor
# return,a
# end
#
#
#
# def convolve(image, psf): #yea i don't get this one...
#/home/vaughan/bitten/SPIRE/smap_pipeline/astrolib/pro/convolve.pro
#     sim = image.shape
#     sc = floor((sim-1)/2)
#     npix = image.size
#
#     conv = npix * np.real()
