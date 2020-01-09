import numpy as np
from scipy.ndimage import shift
import sys
sys.path.append('utilities/')
import gaussian as g
from astropy.io import fits
import matplotlib.pyplot as plt
from math import *
from astropy.convolution import convolve_fft as convolve
import time

def clus_image_model(fwhm, ref_map,x ,y, f):
    #create a psf to later convolve with our map
    psf = g.makeGaussian(24, 24, fwhm=fwhm)

    #pad psf to match the size of our map
    #need to do this otherwise convolution won't work. Doesn't effect the convolution.
    pad_psf = g.padGaussian(ref_map, psf)
    #psf is padded to size of our map so we can define x_size and y_size with that shape.
    y_size, x_size = pad_psf.shape
    #find sources that are inside our map
    bool_arr = [True for i in range(len(f)) if x[i] >= 0 and x[i] <= x_size and y[i] > 0  and y[i] <= y_size]

    #and cut out the sources that are outside our map.
    x = x[bool_arr]
    y = y[bool_arr]
    f = f[bool_arr]

    #create a list of integer values for x, y and find difference between those and the floats
    ix = [ceil(xi) for xi in x]
    dx = [ix[i] - x[i] for i in range(len(x))]
    iy = [ceil(yi) for yi in y]
    dy = [iy[i] - y[i] for i in range(len(y))]

    #image of only sources
    im = np.zeros((y_size, x_size))
    im_s = np.zeros((y_size, x_size))
    for i in range(len(f)):
        s_image = np.zeros((y_size, x_size))
        s_image[iy[i]-1, ix[i]-1] = f[i]
        p_image = shift(s_image, np.array([dx[i],dy[i]]))
        im += p_image


    image = convolve(im, pad_psf)

    return image



if __name__ == '__main__':
    f = np.random.rand(70000)
    x = np.random.rand(70000) * 50
    y = np.random.rand(70000) * 30
    ref_map = np.zeros((30, 50))
    start = time.time()
    print(start, 'start time')
    im = clus_image_model(3, ref_map, x, y, f)
    end = time.time()
    print(end, 'end time')
    print(end - start, 'time elapsed')
    plt.imshow(im)
    plt.show()
