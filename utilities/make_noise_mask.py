################################################################################
# NAME : make_noise_mask
# DATE STARTED : June 11, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE : A support routine for clus_compute_rings that makes a mask if noise is above a threshold.
#
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#           nsims (number of simulations/files to read in)
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################

def make_noise_mask(maps, col):
    """
    Inputs: Maps is the maps object that we get from get_data
            Col is either 0, 1, 2 to specify which band you want.
            0 = PSW, 1 = PMW, 2 = PLW
    """

    for i in range(len(maps)) :
        mapsize = maps[col]['astr'].naxis
        errmap = map[col]['error']
        for j in range(mapsize[0]):
            for k in range(mapsize[1]):
                if maps[col]['error'][j,k] == 0:
                    errmap[j,k] = 0.005
        #call to make a gaussian function and convol it with the map to smooth it out, but I don't think we need that?
                    if map[j,k] >= 0.004:
                        maps[col]['mask'][j,k] = 1
    return maps[col]['mask']
