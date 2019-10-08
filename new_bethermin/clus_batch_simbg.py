
################################################################################
# NAME : clus_get_batch_simbg.py
# DATE STARTED : September 27, 2019
# AUTHORS : Dale Mercado
# PURPOSE : The wrapper for Conley's Bethermin simmulation generator.
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
#
################################################################################

import numpy as np
import sys
sys.path.append('../utilities')
import config

def clus_batch_simbg():

    clusters = ['a0370','a1689','a1835','a2218','a2219','a2390',
              'cl0024','ms0451','ms1054','ms1358','rxj0152','rxj1347']

    clusters = 'a1689' #Is this how you select the one that you want to be working on???

    nclust = len(clusters)

    nsim = 300
    fluxcut == [[10**(-np.arrage(6)+2))],0]

    for iclust in range(nclust):

        # Line for the output save file path

        for isim in range(nsim):

            print('On sim ' + str(isim+1) + 'of' + str(nsims))

            if icut == 0:
                genbethermin == 1
            else:
                genbethermin == 0

            CLUS_SIM_BACKGROUND()


    if saveplot:
        # lines for saving the output files
