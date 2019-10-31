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

def histogram(nsims):

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
                file = 'sz/sim/szout__' + clusname + str(i) + '.json'
                with open(file) as data_file:
                    data = json.load(data_file)
                # data = scipy.io.readsav(config.CLUSDATA + file,python_dict = True)
                # arr_data = list(data.values())[0][k][0]
                Isim[n] = data[k]['increment']
                n += 1
                #this will be changed to:
                #outfile = config.CLUSDATA + outdir + '/sim/' + outname + '_' + params['clusname'][0] + '.json'
                #or
                #outfile = config.CLUSDATA + outdir + '/' + outname + '_' + params['clusname'][0] + '.json'
            histogram_test(10, k, Isim, xbin, clusname, band)
        exit()

def histogram_test(binsize, k, Isim, xbin, clusname, band):
    dI = abs((np.max(Isim) + np.min(Isim)) / binsize)
    # for b in range(binsize):
    #     xbin[b] = np.min(Isim) + (b + 0.5) * dI
    bins = np.linspace(np.min(Isim), np.max(Isim), num=binsize)
    average = np.mean(Isim)
    avg = [average, average]
    with open('sz/sim/szout__' + clusname + 'real.json') as data_file:
        data = json.load(data_file)
    Ireal = data[0]['increment']
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
    plt.show()




if __name__ == '__main__':
    nsims =  50
    histogram(nsims)
