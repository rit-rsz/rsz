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
        if fixed_psf:
            psf_xy = image_shift( psf, x[i] - ix, y[i] -iy, data)
            #call to image shift so need to port that.

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

    # i do not understand the returns on image_shift.

def extend_array(arr, s0, s1, offset): #i literally don't understand what I'm reading here...
    array = make_array(s0, s1, type=)


def subs_to_coord(subs, n_column):
    s = round(subs)
    n = round(n_columns)
    x = s % n
    y = s / n
    return x, y
