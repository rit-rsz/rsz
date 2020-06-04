# python
import numpy as np
import matplotlib.pyplot as plt
from math import *
import sys, time
from astropy.io import fits
from scipy.stats import norm
sys.path.append('../utilities')
from clus_get_data import *
from config import *

def clus_lens_template(name=None,band=None):
    maps,err = clus_get_data('rxj1347','0',sgen=3)

    for i in range(100):
        image_psw = config.OUTPUT + 'pcat_residuals/rxj1347_resid_PSW_%s.fits' %(i)
        image_pmw = config.OUTPUT + 'pcat_residuals/rxj1347_resid_PMW_%s.fits' %(i)
        image_plw = config.OUTPUT + 'pcat_residuals/rxj1347_resid_PLW_%s.fits' %(i)

        # image_psw = 'lense_resid/rxj1347_resid_PSW_%slense_noiseless.fits' %(i)
        # image_pmw = 'lense_resid/rxj1347_resid_PMW_%slense_noiseless.fits' %(i)
        # image_plw = 'lense_resid/rxj1347_resid_PLW_%slense_noiseless.fits' %(i) #%(str(i).zfill(2))
        # image_psw = 'lense_resid/rxj1347_resid_PSW_%slense.fits' %(i)
        # image_pmw = 'lense_resid/rxj1347_resid_PMW_%slense.fits' %(i)
        # image_plw = 'lense_resid/rxj1347_resid_PLW_%slense.fits' %(i) #%(str(i).zfill(2))
        # data = fits.getdata(image_psw)[125,125]
        # bin_cent = create_bin_edges(data, .5)
        # hist,bin = np.histogram(data,bins=bin_cent)
        # mu, sigma = norm.fit(data)
        # y = norm(loc=mu, scale=sigma).pdf(bin_cent)
        # plt.plot(bin_cent, y)
        # plt.hist(data,bin,histtype='step',normed=1)
        # plt.xlabel('mJy/beam')
        # plt.ylabel('# pixels')
        # plt.savefig('lense_hist_test.png')
        # plt.clf()
        # exit()

        image_ob_psw = fits.getdata(image_psw) * maps[0]['mask']
        image_ob_pmw = fits.getdata(image_pmw) * maps[1]['mask']
        image_ob_plw = fits.getdata(image_plw) * maps[2]['mask']

        if i == 0 :
            size_psw = image_ob_psw.shape
            image_list_psw = np.zeros((size_psw[0]*size_psw[1],100))
            image_list_psw[:,0] = image_ob_psw.flatten()

            size_pmw = image_ob_pmw.shape
            image_list_pmw = np.zeros((size_pmw[0]*size_pmw[1],100))
            image_list_pmw[:,0] = image_ob_pmw.flatten()

            size_plw = image_ob_plw.shape
            image_list_plw = np.zeros((size_plw[0]*size_plw[1],100))
            image_list_plw[:,0] = image_ob_plw.flatten()

        else :
            image_list_psw[:,i] = image_ob_psw.flatten()
            image_list_pmw[:,i] = image_ob_pmw.flatten()
            image_list_plw[:,i] = image_ob_plw.flatten()

    map_chunk(name,'PSW',size_psw,image_list_psw)
    map_chunk(name,'PMW',size_pmw,image_list_pmw)
    map_chunk(name,'PLW',size_plw,image_list_plw)

    return

def map_chunk(name,band,size,image_list):

    count = 0
    step = floor((size[0] * size[1]) / 100)
    print('step',step)
    for j in range(step,size[0]*size[1], step):
        if j == step:
            temp = image_list[0:j-1,:]
        else:
            temp = image_list[j-step - 1:j,:]
        np.save(config.HOME + 'Lensing/lense_chunks/lense_chunk_%s_%s_%s.npy'%(name,band,count),temp)
        count += 1

    if (size[0] * size[1] * step) != len(image_list[:,0]):
        temp = image_list[j-1:len(image_list[:,0])-1, :]
        np.save(config.HOME + 'Lensing/lense_chunks/lense_chunk_%s_%s_%s.npy'%(name,band,count),temp)

    print('count:',count)

    return

if __name__ == '__main__' :
    a = time.time()
    clus_lens_template(name='rxj1347',band='PSW')
    b = time.time()
    print('Total Elapsed Time :', b-a)
