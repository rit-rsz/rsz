################################################################################
# NAME : histogram_revisited
# DATE STARTED : May 30, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE : Adapted version of Nick's histogram script in IDL and
# Victoria's / Dale's script in python to create histograms for the probability
# of our output temperature and dI.
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import numpy as np
import scipy.io
import sys
sys.path.append('utilities/')
import config

def temperature_histogram(nsims):

    clusters = ['a0370','a1689','a1835','a2219',
                'cl0024', 'ms0451','ms1054',
                'rxj0152','rxj1347'] # list of clusters we care about
    bands = ['PSW', 'PMW', 'PLW'] #the three different bands

    for j in range(len(clusters)):
        clusname = clusters[j]

        for k in range(len(bands)):
            Isim = np.zeros(nsims)
            xbin = np.zeros(10)
            band = bands[k]

            n = 0

            for i in range(200, 200+nsims):
                file = 'sz/sim/szout_0' + str(i) + '_' + clusname + '.sav'
                #this will be changed to:
                #outfile = config.CLUSDATA + outdir + '/sim/' + outname + '_' + params['clusname'][0] + '.json'
                #or
                #outfile = config.CLUSDATA + outdir + '/' + outname + '_' + params['clusname'][0] + '.json'


if __name__ == '__main__':
    nsims =  1
    temperature_histogram(nsims)
