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
import json
import config
import numpy as np

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

def save_fitsz(increment, offset, radave, params, simflag=1, outname='szout', outdir='sz', verbose=1):
    ncols = len(radave)


    #szout = REPLICATE(szout, ncols)
    szout = []
    for i in range(ncols):
        szout.append(create_template(params))
        szout[i]['increment'] = increment[i].tolist()
        szout[i]['offset'] = offset[i].tolist()
        szout[i]['band'] = radave[i]['band']
        szout[i]['midbin'] = radave[i]['midbin'].tolist()
        szout[i]['fluxbin'] = radave[i]['fluxbin'].tolist()
        szout[i]['errbin'] = radave[i]['errbin'].tolist()

        #i have a feeling this may  not work ex: szout[i]['radbin'] = radave[i]['radbin'] but not sure

    if not simflag:
        outfile = config.HOME + outdir + '/' + outname + '_' + params['clusname'][0] + '.json'
        print(outfile)
    else:
        outfile = config.HOME + outdir + '/sim/' + outname + '_' + params['clusname'][0] + '.json'
        print(outfile)
    with open(outfile, 'w') as f:
        json.dump(szout, f)
