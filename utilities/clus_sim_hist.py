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

def histogram(clusname, nsims, superplot=0):
    bands = ['PSW', 'PMW', 'PLW']
    g = gauss_sample(5, 3, nsims)
    for k in range(len(bands)):
        Isim = np.zeros(nsims)
        xbin = np.zeros(10)
        band = bands[k]
        n = 0
        for i in range(200, 200+nsims):
            # file = config.CLUSDATA + outdir + '/sim/' + outname + '_' + params['clusname'][0] + '.json'
            # with open(file) as data_file:
                # data = json.load(data_file)
            #### for reading old files to test
            ####### data = scipy.io.readsav(config.CLUSDATA + file,python_dict = True)
            ####### arr_data = list(data.values())[0][k][0]
            # Isim[n] = data[k]['increment']
            Isim[n] = g[n]
            n += 1

        histogram_test(10, k, Isim, xbin, clusname, band, superplot)

def gauss_sample(center, std, nsims):
    g = np.random.normal(center, std, nsims)
    return g

def histogram_test(binsize, k, Isim, xbin, clusname, band, superplot=0):
    dI = abs((np.max(Isim) + np.min(Isim)) / binsize)
    # for b in range(binsize):
    #     xbin[b] = np.min(Isim) + (b + 0.5) * dI
    bins = np.linspace(np.min(Isim), np.max(Isim), num=binsize)
    average = np.mean(Isim)
    avg = [average, average]
    inputI = 3
    DI = avg - inputI

    ############### code to put back in
    # inputI = readfromsomefile
    # with open('sz/sim/szout__' + clusname + 'real.json') as data_file:
    #     data = json.load(data_file)
    # Ireal = data[k]['increment']
    # Ir= [Ireal, Ireal]
    if superplot:
        weights = np.ones_like(Isim)/float(len(Isim))
        n, bins, patches = plt.hist(Isim, bins, histtype='step', weights=weights)
        y = [0, np.max(n)]
        # plt.plot(Ir, y, label='Intensity from Real Data %.2E' % Ireal)
        plt.plot(avg, y, label='Average Intensity %.2E' % average)
        plt.xlabel('I_0 (MJy/sr)')
        plt.ylabel('Probability (percent)')
        plt.title('I_0 for %s band: %s' % (clusname, band))
        plt.legend()
        plt.show()

    return DI




if __name__ == '__main__':
    nsims =  23
    histogram('test', 1000)
