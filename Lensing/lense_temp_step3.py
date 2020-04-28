import numpy as np
import sys, time, subprocess, os
sys.path.append('../utilities')
sys.path.append('../source_handling')
import config
from clus_get_data import clus_get_data
from astropy.io import fits
import matplotlib.pyplot as plt
from math import *

bands = ['PSW','PMW','PLW']
maps, err = clus_get_data('rxj1347', 0, resolution='nr', bolocam=None, verbose=1)

for i in range(3):
    size = maps[i]['signal'].shape
    outmap = np.linspace(0,0,num=size[0]*size[1])
    step = floor((size[0] * size[1]) / 100)
    index_arr = np.arange(step,size[0]*size[1],step)

    for k in range(100):
        j = index_arr[k]
        data_or = np.load('/home/vlb9398/rsz/Lensing/lense_chunks/lense_rxj1347_%s_med_%s.npy'%(bands[i],k))
        data = [x*1e-3 for x in data_or]
        if j == step:
            outmap[0:j-1] = data
        else:
            outmap[j-step-1:j] = data

    if k == 100 :
        outmap[j-1:len(outmap)-1] = data

    outmap = outmap.reshape((size[0],size[1]))
    plt.imshow(outmap,origin=0)
    plt.clim([-0.02,0.02])
    plt.colorbar()
    plt.savefig('lense_template_%s.png'%(bands[i]))
    plt.clf()
    hda = fits.PrimaryHDU(outmap,maps[i]['shead'])
    hda.writeto('lense_template_%s.fits'%(bands[i]),overwrite=True)
