################################################################################
# NAME : save_fitsz
# DATE STARTED : June 24, 2019
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

def save_fitsz(increment, offset, radave, params, simflag=1, outname='szout', outdir='sz', verbose=1):
    ncols = len(radave)

    szout_tempplate = {'increment' : 0,
                       'offset' : 0,
                       'band' : ' ',
                       'radbin' : np.empty(radave[0]['radbin'].shape),
                       'midbin' : np.empty(radave[0]['midbin'].shape),
                       'fluxbin' : np.empty(radave[0]['fluxbin'].shape),
                       'errbin' : np.empty(radave[0]['errbin'].shape),
                       'rc' : params['rc'][0],
                       'beta' : params['beta'][0],
                       'clusname' : params['clusname'][0]}
    #szout = REPLICATE(szout, ncols) I don't think my method of np.tiles will work for a dictionary...
    szout = []
    for i in range(ncols):
        szout.append(szout_template)
        szout[i]['increment'] = increment[i]
        szout[i]['offset'] = offset[i]
        szout[i]['band'] = radave[i]['band']
        szout[i]['radbin'] = radave[i]['radbin']
        szout[i]['midbin'] = radave[i]['midbin']
        szout[i]['fluxbin'] = radave[i]['fluxbin']
        szout[i]['errbin'] = radave[i]['errbin']

        #i have a feeling this may  not work ex: szout[i]['radbin'] = radave[i]['radbin'] but not sure

    if not simflag:
        outfile = config.CLUSDATA + outdir + '/' + outname + '_' + params['clusname'][0] + '.sav'
        print(outfile)
    else:
        outfile = config.CLUSDATA + outdir + '/sim/' + outname + '_' + params['clusname'][0] + '.sav'
        print(outfile)
    #call to save the .sav file... maybe we want to make this .sav files .csv files but not sure
