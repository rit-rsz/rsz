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
sys.path.append('utilities/')
import config
import matplotlib.pyplot as plt
import json
from math import *

def clus_sim_hist(nsims,clusname):

    bands = ['PSW','PMW','PLW']
    sz_avg = []
    for k in range(len(bands)):
        Isim = np.zeros(nsims)
        band = bands[k]

        ''' Turned off for testing'''
        # n = 0
        # for i in range(200, 200+nsims):
        #     file = config.HOME + 'output/szout__' + clusname + str(i) + '.json'
        #     with open(file) as data_file:
        #         data = json.load(data_file)
        #     Isim[n] = data[k]['increment']

            # n += 1

        average = histogram_sim(25, k, Isim, clusname, band,nsims)
        sz_avg.append(average)
    return sz_avg

def histogram_sim(binsize, k, Isim, clusname, band, nsims):
    # open the real SZ amplitude for given cluster
    with open(config.HOME + 'outputs/szout/szout__' + clusname + 'real.json') as data_file:
        data = json.load(data_file)

    ''' for testing... load the one sim we have and make it the average dI '''
    with open(config.HOME + 'outputs/szout/szout__' + clusname + '0' + '.json') as data_file:
        sim_data = json.load(data_file)



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
