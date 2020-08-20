
################################################################################
# NAME : test.py
# DATE STARTED : Aug 12, 2020
# AUTHORS : Benjamin Vaughan
# PURPOSE : This is a collection of utility functions for the lensing script.
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
# OUTPUTS :
#
# REVISION HISTORY :
################################################################################import numpy as np
from astropy.convolution import Gaussian2DKernel

def centroid(arr):
    '''
    Purpose: Calculate the weighted average along x and y axes
    Input - array of values
    Output- array of xbar and ybar
    '''
    s = arr.shape
    total = np.sum(arr)

    liny = np.arange(0, s[0])
    linx = np.arange(0, s[1])

    xweight, yweight = np.meshgrid(linx, liny)

    xcm =  np.sum(arr * xweight) / total
    ycm =  np.sum(arr * yweight) / total

    center = np.array([xcm, ycm], dtype=np.float32)
    return center

def find_regions(img, start_x, start_y, max_val, min_val):
    '''
    This function finds connected pixels within a given range (includes diagonals)
    input:
    img - an array representing an image
    start_x - pixel coordinate to start x position at
    start_y - pixel coordinate to start y position at
    max_val - max value for the range
    min_val - min value for the range
    '''

    x_indexes, y_indexes = np.where(img)
    #create a list of all the x_indexes and y_indexes in the image and then switch 0,0 with start_x and start_y
    x_indexes[0] = start_x
    y_indexes[0] = start_y
    x_indexes[start_x] = 0
    y_indexes[start_y] = 0

    npixels = img.size

    #if the image only has 4 pixels we don't need to do any fancy math since all of the pixels are connected together
    if npixels <= 4:
        good_x, good_y = np.where(np.logical_and(img >= min_val, img <= max_val))


    #otherwise we have to do some ugly math
    else:
        count = 0 #initialize a counter for the while loop
        #start with our initial x and y values
        good_x = [start_x]
        good_y = [start_y]
        good = [[start_x, start_y]]

        #condition that stops once we run out of connected pixels in our given range
        while count <= len(good_x) - 1:

            #initialize subscripts to look at a 3x3 patch centered on the pixel of interest
            min_x = good_x[count] - 1
            max_x = good_x[count] + 2
            min_y = good_y[count] - 1
            max_y = good_y[count] + 2

            #this is a bunch of checks to see if our patch is going over the edge of image.
            if min_x < 0:
                min_x = 0
            if max_x > img.shape[0]:
                max_x = good_x[count] + 1
                if max_x > img.shape[0]:
                    max_x = good_x[count]
            if min_y < 0:
                min_y = 0
            if max_y > img.shape[1]:
                max_y = good_y[count] + 1
                if max_y > img.shape[1]:
                    max_y = good_y[count]

            #pull the patch out of the image
            patch = img[min_x:max_x, min_y:max_y]

            #find values where z is in range of min_val <= z <= max_val
            true_x, true_y = np.where(np.logical_and(patch >= min_val, patch <= max_val))

            #iterate our counter + 1
            count += 1
            for indx, indy in zip(true_x, true_y):
                #add min_x and min_y for our patch back so we get the indexes for the whole image
                new_x_ind = indx + min_x
                new_y_ind = indy + min_y
                #check for double counting
                if [new_x_ind,new_y_ind] not in good:
                    #append our x, y indexes
                    good.append([new_x_ind, new_y_ind])
                    good_x.append(new_x_ind)
                    good_y.append(new_y_ind)

            #back up check
            #if we have iterated over every pixel then something is wrong.
            if count > npixels:
                break

        #this is to check that our start_x and start_y have z in range of min_val <= z <= max_val
        first_val = img[good_x[0], good_y[0]]
        if not np.logical_and(first_val >= min_val, first_val <= max_val):
            good_x.remove(good_x[0])
            good_y.remove(good_y[0])

    return good_x, good_y


def deriv(y):
    #This ia python port of IDLS DERIV function found at https://www.harrisgeospatial.com/docs/DERIV.html
    '''
    This function uses a three-point lagrangian interpolation to compute the derivative of an evenenly-spaced array of data.
    Inputs:
    y - y values of an array
    Outputs:
    d - the derivative of y with respect to x
    '''
    n = len(y)
    d = np.subtract( np.roll(y, -1), np.roll(y, 1)) / 2.
    d[0] = (-3.0 *y[0] + 4.0 * y[1] - y[2]) / 2.
    d[n-1] = (3.0 * y[n-1] - 4.0*y[n-2] + y[n-3]) / 2.
    return d

def get_gauss_beam(fwhm, pixscale, band, nfwhm=5.0, oversamp=1):
    retext = round(fwhm * nfwhm / pixscale)
    if retext % 2 == 0:
        retext += 1

    bmsigma = fwhm / math.sqrt(8 * math.log(2))

    beam = Gaussian2DKernel(bmsigma / pixscale, x_size=retext,
                            y_size=retext, mode='oversample',
                            factor=oversamp)
    beam *= 1.0 / beam.array.max()
    return beam

def get_gaussian_psf_template(pixel_fwhm, nbin=5):
    nc = nbin**2
    psfnew = Gaussian2DKernel(pixel_fwhm/2.355*nbin, x_size=125, y_size=125).array.astype(np.float32)
    psfnew2 = psfnew / np.max(psfnew) * nc
    cf = psf_poly_fit(psfnew, nbin=nbin)
    return psfnew2, cf, nc, nbin


import numpy as np


def psf_poly_fit(psf0, nbin):
        assert psf0.shape[0] == psf0.shape[1] # assert PSF is square
        npix = psf0.shape[0]

        # pad by one row and one column
        psf = np.zeros((npix+1, npix+1), dtype=np.float32)
        psf[0:npix, 0:npix] = psf0

        # make design matrix for each nbin x nbin region
        # print(type(npix), type(nbin))
        nc = int(npix/nbin) # dimension of original psf
        nx = nbin+1
        y, x = np.mgrid[0:nx, 0:nx] / np.float32(nbin)
        x = x.flatten()
        y = y.flatten()
        A = np.column_stack([np.full(nx*nx, 1, dtype=np.float32), x, y, x*x, x*y, y*y, x*x*x, x*x*y, x*y*y, y*y*y]).astype(np.float32)
        # output array of coefficients

        cf = np.zeros((A.shape[1], nc, nc), dtype=np.float32)

        # loop over original psf pixels and get fit coefficients
        for iy in range(nc):
            for ix in range(nc):
                # solve p = A cf for cf
                p = psf[iy*nbin:(iy+1)*nbin+1, ix*nbin:(ix+1)*nbin+1].flatten()
                AtAinv = np.linalg.inv(np.dot(A.T, A))
                ans = np.dot(AtAinv, np.dot(A.T, p))
                cf[:,iy,ix] = ans

        return cf.reshape(cf.shape[0], cf.shape[1]*cf.shape[2])

def image_model_eval(x, y, f, back, imsz, nc, cf, regsize=None, margin=0, offsetx=0, offsety=0, weights=None, ref=None, lib=None):

    # assert f.dtype == np.float32
    # not sure what to do with cf
    #assert cf.dtype == np.float32
    if ref is not None:
        assert ref.dtype == np.float32

    if weights is None:
        weights = np.full(imsz, 1., dtype=np.float32)
    # elif len(weights)==1:
    #     weights = np.full(imsz, weights, dtype=np.float32)
    if regsize is None:
        regsize = max(imsz[0], imsz[1])

    # FIXME sometimes phonions are outside image... what is best way to handle?
    goodsrc = (x > 0) * (x < imsz[0] - 1) * (y > 0) * (y < imsz[1] - 1)
    x = x.compress(goodsrc)
    y = y.compress(goodsrc)
    f = f.compress(goodsrc)

    nstar = x.size
    rad = round(nc/2) # 12 for nc = 25

    nregy = int(imsz[1]/regsize + 1) # assumes imsz % regsize = 0?
    nregx = int(imsz[0]/regsize + 1)

    ix = np.ceil(x).astype(np.int32)
    dx = ix - x
    iy = np.ceil(y).astype(np.int32)
    dy = iy - y

    dd = np.column_stack((np.full(nstar, 1., dtype=np.float32), dx, dy, dx*dx, dx*dy, dy*dy, dx*dx*dx, dx*dx*dy, dx*dy*dy, dy*dy*dy)).astype(np.float32) * f[:, None]
    if lib is None:
        image = np.full((imsz[1]+2*rad+1,imsz[0]+2*rad+1), back, dtype=np.float32)
        recon2 = np.dot(dd, cf).reshape((nstar,nc,nc))
        recon = np.zeros((nstar,nc,nc), dtype=np.float32)
        recon[:,:,:] = recon2[:,:,:]
        for i in range(nstar):
            image[iy[i]:iy[i]+rad+rad+1,ix[i]:ix[i]+rad+rad+1] += recon[i,:,:]

        image = image[rad:imsz[1]+rad,rad:imsz[0]+rad]

        if ref is not None:
            diff = ref - image
            diff2 = np.zeros((nregy, nregx), dtype=np.float64)
            for i in range(nregy):
                y0 = max(i*regsize - offsety - margin, 0)
                y1 = min((i+1)*regsize - offsety + margin, imsz[1])
                for j in range(nregx):
                    x0 = max(j*regsize - offsetx - margin, 0)
                    x1 = min((j+1)*regsize - offsetx + margin, imsz[0])
                    subdiff = diff[y0:y1,x0:x1]
                    diff2[i,j] = np.sum(subdiff*subdiff*weights[y0:y1,x0:x1])
    else:
        image = np.full((imsz[1], imsz[0]), back, dtype=np.float32)
        recon = np.zeros((nstar,nc*nc), dtype=np.float32)
        reftemp = ref
        if ref is None:
            reftemp = np.zeros((imsz[1], imsz[0]), dtype=np.float32)
        diff2 = np.zeros((nregy, nregx), dtype=np.float64)
        lib(imsz[0], imsz[1], nstar, nc, cf.shape[0], dd, cf, recon, ix, iy, image, reftemp, weights, diff2, regsize, margin, offsetx, offsety)

    # print image.shape, nc
    if ref is not None:
        return image, diff2
    else:
        return image
