################################################################################
# NAME : save_fitsz
# DATE STARTED : June 24, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE : save sz data to a .json file
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import config
import numpy as np
import os

def save_fitsz(increment, offset, params, sgen=None, nsim=0, verbose=1):
    szout = []

    for i in range(3):
        szout_template = {'increment' : increment[i],
                          'offset' : offset[i],
                          'rc' : params['rc'],
                          'beta' : params['beta'],
                          'clusname' : params['clusname']}
                           # 'band' : radave[i]['band'],
                           # 'midbin' : radave[i]['midbin'],
                           # 'fluxbin' : radave[i]['fluxbin'],
                           # 'errbin' : radave[i]['errbin'],

        szout.append(szout_template)

    if sgen is not None:
        outfile = config.OUTPUT + 'szout/' + params['clusname'] + 'szout_' + str(nsim) + '.npy'
    else:
        outfile = config.OUTPUT + 'szout/' + params['clusname'] + 'szout_' + 'real.npy'
    print('save fitsz has saved the file')
    np.save(outfile,szout,allow_pickle=True)
