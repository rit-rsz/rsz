################################################################################
# NAME : subtract_cat.py
# DATE STARTED : June 21, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE : idk lol
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import numpy as np

def subtract_cat(maps, cat, verbose=1):
    err = False
    ncols = len(maps)

    for i in range(ncols):
        if verbose:
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
            datasub[value] = maps[i]['signal'].data[value] - cat[i]['model'][value]
        datafull = maps[i]['signal'].data

        whpl = []
        for j in range(maps[i]['mask'].shape[0]):
            for k in range(maps[i]['mask'].shape[1]):
                if maps[i]['mask'] == 1:
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
