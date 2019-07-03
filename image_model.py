################################################################################
# NAME : image_model.py
# DATE STARTED : June 27, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE : create an image based off of some data.
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import numpy as np

def image_model(x,y,f,x_size,y_size, psf, data, reference_pix=0,
                extra, lx, ux, ly, uy, model):
    """
    Inputs: x -
            y -
            f - source flux an array of fluxes
            x_size -
            y_size -
            psf - the psf for the band we are looking at.
    """
    nstar = len(f) # f is the flux, nstar is number of stars

    fixed_psf = len(psf.shape)

    psf_size = psf.shape
    if not reference_pix:
        psf_ref_pix = get_max(psf[:,:,0]) #i think this is the python equivlanet of psf[*,*,0]

    for i in range(nstar):
        ix = round(x[n])
        iy = round(y[n])
        psf_xy = image_shift( psf, x[i] - ix, y[i] -iy, data)
        image = gaussian_psf()

def gaussian_psf()

def get_max(array):
    s = array.shape
    x = False
    if len(s) != 2:
        return False
    #the way that they coded starfinder is kind of hard to follow, but i think this is correct...
    indexes = []
    if all:
        for i in range(s[0]):
            for j in range(s[1]):
                if s[i,j] >= np.amax(array):
                    indexes.append([i,j])
    else:
        m = np.amax(array)
        indexes = np.where(array == m)
    x, y = subs_to_coord(indexes, s[0])
    p = [x[0]], y[0]]
    return p

def image_shift(image, x_shift, y_shift, data, interp_type=None):
    size = image.shape
    extsize = size + 2
    imag, offset = extend_array(image, extsize[0], extsize[1])
    lo = offset
    up = lo + size - 1

    return imag[lo[0]:up[0], lo[1]:up[1]]

def extend_array(arr, s0, s1, offset):
    s = array.shape
    if s0 == s[0] and s1 == s[1] or s0 =< s[0] or s1 =< s[1]:
        array1 = array
        o = [0,0]
    else:
        o = [s0,s1] - s
        o = (o + o % 2) / 2
    array1 = np.empty(s0, s1)
    array1[o[0], o[1]] = arr
    return array1, o


def subs_to_coord(subs, n_column):
    s = round(subs)
    n = round(n_columns)
    x = s % n
    y = s / n
    return x, y
