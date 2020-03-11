################################################################################
# NAME : histogram_revisited
# DATE STARTED : October 23, 2019
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
sys.path.append('../utilities/')
import config
import matplotlib.pyplot as plt
import json
from math import *

def clus_sim_hist(nsims,clusname):
    tot_sims = 10
    bands = ['PSW','PMW','PLW']
    sz_avg = []
    for k in range(len(bands)):
        Isim = np.zeros(tot_sims)
        band = bands[k]
        average = histogram_sim(25, k, Isim, clusname, band, nsims)
        sz_avg.append(average)
    return sz_avg

def histogram_sim(binsize, k, Isim, clusname, band, nsims):
    # open the real SZ amplitude for given cluster
    data_file = config.OUTPUT + 'szout/' + clusname + 'szout_' + 'real.npy'
    data = np.load(data_file)

    data_file_sim = config.OUTPUT + 'szout/' + clusname + 'szout_' + str(nsims) + '.npy'
    sim_data = np.load(data_file_sim)

    Ireal = data[k]['increment']
    Ifake = sim_data[k]['increment']
    Isim = np.random.normal(loc=Ifake,scale=3.0 ,size=nsims)

    dI = abs((np.max(Isim) + np.min(Isim)) / binsize)
    bins = np.linspace(np.min(Isim), np.max(Isim), num=binsize)
    average = np.mean(Isim)
    avg = [average, average]
    Ir= [Ireal, Ireal]
    weights = np.ones_like(Isim)/float(len(Isim))
    n, bins, patches = plt.hist(Isim, bins, histtype='step', weights=weights)
    y = [0, np.max(n)]
    plt.plot(Ir, y, label='Intensity from Real Data %.2E' % Ireal)
    plt.plot(avg, y, label='Average Intensity %.2E' % average)
    plt.xlabel('I_0 (MJy/sr)')
    plt.ylabel('Probability (percent)')
    plt.title('I_0 for %s band: %s' % (clusname, band))
    plt.legend()
    plt.savefig(config.HOME + 'outputs/dI_hist/dI_hist_%s_%s_TEST.png' %(clusname,band))
    plt.clf()

    return average

if __name__ == '__main__':
    nsims =  99
    clusname = 'rxj1347' # list of clusters we care about
    clus_sim_hist(nsims,clusname)
