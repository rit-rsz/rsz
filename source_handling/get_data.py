
################################################################################
# NAME : get_data.py
# DATE STARTED : June 11, 2019
# AUTHORS : Victoria Butler & Dale Mercado
# PURPOSE : This is where the data is retevied and can be
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import scipy.io
import numpy as np
# from config import * #(this line will give me access to all directory variables)
import matplotlib.pyplot as plt
import math
from astropy.io import fits
import os

def get_data(clusname, verbose = 1):
    verbose = []
    resolution = []
    version = []
    bolocam = []
    manpath = 0

    if not verbose:
        verbose = 1
    else:
        verbose = verbose

    if not resolution:
        resolution = 'nr'
    else:
        resolution = resolution

    if not version:
        version = '1'
    else:
        version = version

    if not bolocam:
        cols = ['PSW','PMW','PLW']
        bolocam = 0
    else:
        cols = ['PSW','PMW','PLW','BOLOCAM']
        bolocam = 1

    if manpath == 0:

        hermfile = [x for x in x.startswith('/data/mercado/SPIRE/' + 'hermes_cluster/' + clusname)\
                    + x.endswith('_' + resolution + '_' + version + '.fits')]
        # This is still broken and I am usure how to move forward on this
        print(hermfile)
        # hermfile = str('/data/mercado/SPIRE/' + 'hermes_cluster/' + clusname + '*_' +\
                    # resolution + '_' + version + '.fits')

        hermfiles = fits.open(hermfile)
        print(hermfiles[0])

        hlsfile = (CLUSDATA + 'hls_cluster/' + clusname +'.fits')
        hlsfiles = fits.open()

        snapfile = (CLUSDATA + 'snapshot_cluster/' + clusname +'*.fits')
        snapfiles = fits.open()

        if snapcount ^ hermcount ^ hlscount != 3 :
            errmsg = ('Problem finding exclusive files, file dump is:', \
                        snapfiles, hermfiles, hlsfiles)
            if verbose:
                print errmsg
            return success
        if hermcount == 3:
            files = hermfiles
            nfiles = hermcount

        if hlscount == 3:
            files = hlsfiles
            nfiles = hlscount

        if snapcount == 3:
            files = snapfiles
            nfiles = snapcount
    else:
        files # Need to apply the file search again
        if mancount == 3:
            errmsg = 'Cannot find 3 files fitting that path description!'
            #message ERRMSG
            return success
        nfiles = mancount

    if bolocam:
        nfiles = nfiles +1

#     for i in range(nfiles, nfiles-1)
# #       Need to get this to look more like the idl code
# #       There it runs through ifile= 0 to nfiles -1
#         if i < 3:
#             errmsg



# def read_file():
#
#     calfac = some Number
#     # This will ultimatly be in the list of constants
#
#
#     # The rest of the scrpit involves idl_libs stuff that
#     # will get grabbed from astropy






if __name__ == '__main__':
    get_data('rxj1347',1)
