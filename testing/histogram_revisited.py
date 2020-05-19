################################################################################
# NAME : histogram_revisited
# DATE STARTED : October 23, 2019
# AUTHORS : Benjamin Vaughan & Victoria Butler
# PURPOSE : Adapted version of Nick's histogram script in IDL and
# Victoria's / Dale's script in python to create histograms for the probability
# of our output temperature and dI.
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
# OUTPUTS :
# REVISION HISTORY : VLB - 04/29/20 : changed to using numpy szout files
################################################################################
import numpy as np
import sys
sys.path.append('../utilities')
import config
import matplotlib.pyplot as plt
from math import *

def histogram(nsims,cluster):

    bands = ['PSW', 'PMW', 'PLW'] #the three different bands

    for k in range(len(bands)):
        Isim = []
        band = bands[k]

        for i in range(nsims):
            file = config.OUTPUT + 'szout/' + str(cluster) + 'szout_' + str(i) + '.npy'
            data = np.load(file,allow_pickle=True)
            Isim.append(data[k].get('increment'))

        histogram_test(50, k, Isim, cluster, band)
    return

def histogram_test(binsize, k, Isim, clusname, band):
    dI = abs((np.max(Isim) + np.min(Isim)) / binsize)
    bins = np.linspace(np.min(Isim), np.max(Isim), num=binsize)
    average = np.mean(Isim)
    avg = [average, average]
    data = np.load(config.OUTPUT + 'szout/' + str(clusname) + 'szout_real.npy',allow_pickle=True)
    Ireal = data[k].get('increment')
    Ir = [Ireal, Ireal]

    weights = np.ones_like(Isim)/float(len(Isim))
    n, bins, patches = plt.hist(Isim, bins, histtype='step', weights=weights)
    y = [0, np.max(n)]
    plt.plot(Ir, y, label='Real Data dI %.2E' % Ireal)
    plt.plot(avg, y, label='Average dI %.2E' % average)
    plt.xlim(-1.5,1.5)
    plt.xlabel('I_0 (MJy/sr)')
    plt.ylabel('Probability (percent)')
    plt.title('I_0 for %s band: %s' % (clusname, band))
    plt.legend(loc='upper right')
    plt.savefig(config.OUTPUT + 'dI_hist/dI_hist_%s_%s.png'%(clusname,band))
    plt.clf()


if __name__ == '__main__':
    # clusters = ['a0370','a1689','a1835','a2219','cl0024', 'ms0451','ms1054','rxj0152','rxj1347']
    histogram(100,'rxj1347')
