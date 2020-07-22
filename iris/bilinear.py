################################################################################
# Name : mosaiq
# Purpose : This does a transformation on the input image through a
# bilinear interpolation method
#Author : Benjamin Vaughan
#Start Date : Oct 4, 2019
#Additional Info
#
################################################################################
import numpy as np
from scipy.interpolate import interpn as interp
import matplotlib.pyplot as plt
from utility import bilinear_interpolation

def mbilinear(x, y, array):
    array = np.squeeze(array)
    six = x.shape
    siy = y.shape
    sia = array.shape
    Nax = sia[0]
    Nay = sia[1]
    Nx = six[0]
    Ny = six[1]
    output = np.zeros((Nx, Ny))
    minval = np.min(array)


    count = 0
    for i in range(Nx):
        for j in range(Ny):
            if x[i,j] < 0 or x[i,j] > Nax-1 or y[i,j] < 0 or y[i,j] > Nay-1:
                count += 1

    inter_percent = 1. * (Nx * Ny - count) / (Nx * Ny) * 100
    print('Images Intersection = %s percent' % (inter_percent))

    pixx = np.arange(0, Nax)
    pixy = np.arange(0, Nay)

    for i in range(Ny):
        mind = []
        for j in range(Nx):
            if x[j,i] >= 0 and x[j,i] <= Nax-1 and y[j,i] >= 0 and y[j,i] <= Nay-1:
                mind.append(j)
        mind = np.asarray(mind)
        if len(mind) != 0:
            xx = x[mind, i]
            yy = y[mind, i]
            trunc = bilinear_interpolation(xx, yy, array)
            output[mind, i] = trunc[:]

    ind = []
    for i in range(output.shape[0]):
        for j in range(output.shape[1]):
            if minval > output[i,j]:
                output[i,j] = -32768

    return output
