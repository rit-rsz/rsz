################################################################################
# NAME : subtract_cat.py
# DATE STARTED : June 21, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE : the purpose of this script is to subtract the data from XID from
#           the original signal
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import numpy as np
import matplotlib.pyplot as plt
from get_data import *

def subtract_cat(maps, cat, verbose=1):
    err = False
    ncols = len(maps)
    print(len(maps), 'length of initial maps')
    for i in range(ncols):
        if verbose:
            print(len(maps[i]))
            print('Subtracting for %s' %(maps[i]['band']))
        datafull = np.empty(maps[i]['astr']['NAXIS'])
        datasub = np.empty(maps[i]['astr']['NAXIS'])
        counter = 0
        whgd = []
        # whpl = np.where(np.isfinite(maps[i]['signal'] == False))
        for j in range(maps[i]['signal'].data.shape[0]):
            for k in range(maps[i]['signal'].data.shape[1]):
                if np.isfinite(maps[i]['signal'].data[j,k]) == False:
                    pass
                else:
                    whgd.append([j,k])
        for value in whgd:
            datasub[value] = maps[i]['signal'].data[value] - cat[i]['signal'].data[value] # we used cat[i]['signal'].data to test but we should really be using cat[i]['model'] in the actual version.
        datafull = maps[i]['signal'].data

        # print(np.max(datasub), np.min(datasub))
        # plt.imshow(datasub)
        # plt.show()

        whpl = []
        for j in range(maps[i]['mask'].shape[0]):
            for k in range(maps[i]['mask'].shape[1]):
                if maps[i]['mask'][j,k] == 1:
                    whpl.append([j,k])
        # whpl = np.where(maps[i]['mask'] == 1)
        for value in whpl:
            datasub[value] = np.nan

     # CONTOUR,datafull,/NODATA,$
     #         TITLE='clus_subtract_cat: Signal map for ' + $
     #         (*maps[icol]).band
     # tvimage,bytscl(datafull,min=-0.01,max=0.01,/NAN),/OVERPLOT
     #
     # CONTOUR,datasub,/NODATA,$
     #         TITLE='clus_subtract_cat: Catalog subtracted map for ' + $
     #         (*maps[icol]).band
     # tvimage,bytscl(datasub,min=-0.01,max=0.01,/NAN),/OVERPLOT

        maps[i]['scrrm'] = datasub
        maps[i]['xclean'] = datasub

    return maps, err

if __name__ == "__main__":
    maps, err = get_data('rxj1347')
    maps = subtract_cat(maps, maps)
