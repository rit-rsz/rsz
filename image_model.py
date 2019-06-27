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

    if fixed_psf != 2:
        space_var = psf.shape
        #not really sure what the rest of the code is doing after this

    #this is weird because we can just call ndarray.shape but whatever.
    if fixed_psf or space_var:
        psf_size = psf.shape
        if not reference_pix:
            psf_ref_pix = get_max(psf[:,:,0]) #i think this is the python equivlanet of psf[*,*,0]

    for i in range(nstar):
        ix = round(x[n])
        iy = round(y[n])
        if fixed_psf:
            psf_xy = #call to image shift so need to port that.

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








def subs_to_coord(subs, n_column):
    s = round(subs)
    n = round(n_columns)
    x = s % n
    y = s / n
    return x, y
