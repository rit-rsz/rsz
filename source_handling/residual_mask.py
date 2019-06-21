################################################################################
# NAME : residual_mask.py
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
def residual_mask(maps, pswonly=1, verbose=1):
    err = False
    if pswonly:
        ncols = 1
    else:
        ncols = 3

    for i in range(ncols):
        snmap = np.empty(maps[i]['srcrm'].shape)
        for j in range(maps[i]['srcrm'].shape[0]):
            for k in range(maps[i]['srcrm'].shape[1]):
                snmap[j,k] = maps[i]['srcrm'][j,k] / maps[i]['error'][j,k]
        abs_snmap = np.empty(snmap.shape) #i think these are 2 dimensional arrays...
        for j in range(snmap.shape[0]):
            for k in range(snmap.shape[1]):
                abs_snmap[j,k] = abs(snmap[j,k])

        #the stuff after this is commented out.
        #then it gives a call to contour...
        maskmap = np.empty(maps[i]['mask'].shape)
        for j in range(maps[i]['mask'].shape[0]):
            for k in range(maps[i]['mask'].shape[1]):
                    maskmap[j,k] = abs(1 - maps[i]['mask'][j,k])

    return maskmap, err
