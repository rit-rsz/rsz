import numpy as np
import sys, time, subprocess, os
sys.path.append('../utilities')
sys.path.append('../source_handling')
import config
from clus_get_data import clus_get_data
from clus_make_noise_mask import clus_make_noise_mask
from astropy.io import fits
import matplotlib.pyplot as plt
import math

bands = ['PSW','PMW','PLW']
maps, err = clus_get_data('rxj1347', 0, resolution='nr', bolocam=None, verbose=1)

for i in range(3):
    size = maps[i]['signal'].shape
    error = maps[i]['error'].flatten()
    mask = maps[i]['mask'].flatten()
    np.random.seed(102793)
    noise = np.random.normal(loc=0.0,scale=1.0,size=(size[0],size[1])).flatten()
    mask = clus_make_noise_mask(maps,i)

    for j in range(len(error)):
        if not math.isnan(error[j]) :
            noise[j] = error[j]*noise[j]

    for k in range(100):
        data = fits.getdata('/home/vlb9398/rsz/Lensing/lense_resid/rxj1347_resid_%s_%slense.fits'%(bands[i],k))
        outmap = data.flatten() - noise
        outmap = outmap.reshape((size[0],size[1]))*mask

        # plt.imshow(outmap,origin=0)
        # plt.clim([-0.005,0.005])
        # plt.colorbar()
        # plt.savefig('noiseless_lense_template_%s.png'%(bands[i]))
        # plt.clf()
        # exit()

        hda = fits.PrimaryHDU(outmap,maps[i]['shead'])
        hda.writeto('/home/vlb9398/rsz/Lensing/lense_resid/rxj1347_resid_%s_%slense_noiseless.fits'%(bands[i],k),overwrite=True)
