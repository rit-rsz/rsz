import argparse
import numpy as np
import matplotlib.pyplot as plt
from math import *
import sys, os
import time
from scipy.stats import norm

def make_pix_hist(name,chunk_num,band):

    image_list = np.load('/home/vlb9398/rsz/Lensing/lense_chunks/lense_chunk_%s_%s_%s.npy'%(name,band,chunk_num))
    bin_size = 1e-1
    size = np.array(image_list).shape[0]
    mean_image = []
    image = np.nan_to_num(image_list)
    time1 = time.time()
    for j in range(size):
        data = [x*1e3 for x in image[j,:]]
        if sum(data) == 0.0 :
            print('all zeros')
        else :
            bin_cent = create_bin_edges(data, .5)
            hist,bin = np.histogram(data,bins=bin_cent)
            mu, sigma = norm.fit(data)
            y = norm(loc=mu, scale=sigma).pdf(bin_cent)
            # bin_middles = (bin[:-1] + bin[1:]) / 2.0
            mean_image.append(mu)

            # plt.plot(bin_cent, y)
            # plt.hist(data,bin,histtype='step',normed=1)
            # plt.savefig('lense_hist_test.png')
            # plt.clf()

    np.save('/home/vlb9398/rsz/Lensing/lense_chunks/lense_%s_%s_med_%s.npy'%(name,band,chunk_num),mean_image)
    os.remove('/home/vlb9398/rsz/Lensing/lense_chunks/lense_chunk_%s_%s_%s.npy'%(name,band,chunk_num))
    time2 = time.time()
    print('Total sub-map runtime %s:'%(chunk_num),time2-time1)

    return

def create_bin_edges(pset, bin_size):
    #the purpose of this function is to create the binning for our pixels
    '''
    pset = the set of pixels at a given coordinate in each image
    bin_size = bin_size for the histogram
    '''
    min_pix = np.min(pset) #find max and min pixel values for limits of histogram
    max_pix = np.max(pset)
    N = ceil((max_pix - min_pix) / bin_size)
    stop = min_pix + (bin_size * N)
    bin_cent, stepcheck = np.linspace(min_pix, stop, N+1, retstep=True) #create bins based on the uncertainty in pixel vals
    return bin_cent

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-run", help="Process Median of Lensing Map Chunk Pixels", nargs=3, metavar=('name','chunk_num','band'))
    args = parser.parse_args()
    if args.run:
        print(args.run)
        name = args.run[0]
        chunk_num = int(args.run[1]) - 1
        band = args.run[2]
        make_pix_hist(name,chunk_num,band)
