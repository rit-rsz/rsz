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

def create_template(params):
    szout_template = {'increment' : 0,
                       'offset' : 0,
                       'band' : ' ',
                       'midbin' : None,
                       'fluxbin' : None,
                       'errbin' : None,
                       'rc' : params['rc'],
                       'beta' : params['beta'],
                       'clusname' : params['clusname']}
    return szout_template

def save_fitsz(increment, offset, radave, params, sgen=None, nsim=0, verbose=1):
    ncols = len(radave)
    szout = []

    for i in range(ncols):
        szout.append(create_template(params))
        szout[i]['increment'] = increment[i]
        szout[i]['offset'] = offset[i]
        szout[i]['band'] = radave[i]['band']
        szout[i]['midbin'] = radave[i]['midbin']
        szout[i]['fluxbin'] = radave[i]['fluxbin']
        szout[i]['errbin'] = radave[i]['errbin']

    if sgen is not None:
        outfile = config.OUTPUT + 'szout/' + params['clusname'] + 'szout_' + str(nsim) + '.npy'
    else:
        outfile = config.OUTPUT + 'szout/' + params['clusname'] + 'szout_' + 'real.npy'

    np.save(outfile,szout,allow_pickle=True)
