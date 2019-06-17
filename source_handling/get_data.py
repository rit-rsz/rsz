
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
    # place holders for the if statements to work until I add them to the input for get data
    # This will happen when the script is fully functional
    verbose = []
    resolution = []
    version = []
    bolocam = []
    manpath = 0
    CLUSDATA = '/data/mercado/SPIRE/'

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
        dir = []
        # hermfile = [self.dir + x for x in x.startswith('/data/mercado/SPIRE/' + 'hermes_cluster/' + clusname)\
        #             + x.endswith('_' + resolution + '_' + version + '.fits')]
        # files = [self.dir + x for x in os.listdir(self.dir) if x.startswith("omnilog")]
        # This is still broken and I am usure how to move forward on this
        # print(hermfile)
        # hermfile = str('/data/mercado/SPIRE/' + 'hermes_clusters/' + clusname + '*_' +\
        #             resolution + '_' + version + '.fits')
        # This is an example of the file name that we are trying to call
        # rxj1347_PLW_fr_1.fits
        # hermfiles = []
        for i in range(len(cols)):
            # Note Need to change all the hard coded
            hermfile = (CLUSDATA + 'hermes_clusters/' + clusname + '_' +\
                        cols[i] + '_' + resolution + '_' + version + '.fits')
            hermfiles = (fits.open(hermfile))
            # In the idl file it uses a * to grab all the bands here I had to use a for loop in order
            # to acheive the same results
            # This method provides difficulty for the hls and snapfiles below as they use other numbers in their names

            hlsfile = (CLUSDATA + 'hls_clusters/' + clusname + '*')
            hlsfiles = fits.open(hlsfile)
            print(hlsfiles.info())

            snapfile = (CLUSDATA + 'snapshot_clusters/' + clusname +'*.fits')
            snapfiles = fits.open(snapfile)

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



def read_file():
#
#
    calfac = (np.pi/180) * (1/3600)**2 * (np.pi / (4 * np.log(2))) * (1e6)

#
# This will ultimatly be in the list of constants
#
#
#     # The rest of the scrpit involves idl_libs stuff that
#     # will get grabbed from astropy






if __name__ == '__main__':
    get_data('rxj1347',1)
