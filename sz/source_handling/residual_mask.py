################################################################################
# NAME : residual_mask.py
# DATE STARTED : June 21, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE : this function collects the residual mask from maps input in.
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import matplotlib as plt
import matplotlib.pyplot as p
import sys
sys.path.append('source_handling')
sys.path.append('../utilities')
from get_data import *
import numpy as np
from writefits import *
import matplotlib.pyplot as plt

def residual_mask(maps, pswonly=0, verbose=1):
    err = False
    if pswonly:
        ncols = 1
    else:
        ncols = 3


    for i in range(ncols):
        #generate signal to noise maps
        srcrm = maps[i]['srcrm']
        error = maps[i]['error'].data
        snmap = np.empty(srcrm.shape)
        for j in range(srcrm.shape[0]):
            for k in range(srcrm.shape[1]):
                snmap[j,k] = srcrm[j,k] / error[j,k]
        #find s/n pixels > 8
        # whpl = []
        # for j in range(snmap.shape[0]):
        #     for k in range(snmap.shape[1]):
        #         if abs(snmap[j,k]) >= 3:
        #             whpl.append([j,k])



        #the stuff after this is commented out.

        #then it gives a call to contour...
        #and finally make an image I can look at

        plt.imshow(snmap)
        plt.show()



        # contoured = plt.contour(snmap) #maybe this can be a replacement? what does tvimage do in the idl version.
        # p.clabel(contoured, inline=True)
        # p.title('clus_residual_mask : Residual bright source subtracted map for %s' % (maps[i]['band']))
        # p.show() # for some reason this only shows the first band.
        maskmap = np.empty(maps[i]['mask'].shape)
        data = np.empty(maskmap.shape)
        for j in range(maps[i]['mask'].shape[0]):
            for k in range(maps[i]['mask'].shape[1]):
                    maskmap[j,k] = abs(1 - maps[i]['mask'][j,k])
                    data[j,k] = maskmap[j,k] * maps[i]['xclean'][j,k]

        writefits('test_' + str(maps[i]['band']) + '.fits', data=data, header_dict=maps[i]['shead'])

    return maskmap, err


if __name__ == '__main__':
    psw = fits.open('../fits_files/xid_9_subtracted_a0370_PSW.fits')
    pmw = fits.open('../fits_files/xid_9_subtracted_a0370_PMW.fits')
    plw = fits.open('../fits_files/xid_9_subtracted_a0370_PLW.fits')
    maps, err = get_data('a0370')
    maps[0]['srcrm'] = psw[0].data
    maps[1]['srcrm'] = pmw[0].data
    maps[2]['srcrm'] = plw[0].data
    mask, err = residual_mask(maps)
