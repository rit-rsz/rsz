
################################################################################
# NAME : histogram_test.py
# DATE STARTED : May 30, 2019
# AUTHORS : Victoria Butler & Dale Mercado
# PURPOSE : Create a python couterpart to the script sz_KS_nick_3.pro.
#           This script will be used to test the creation of Nick's Histograms
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import scipy.io
import numpy as np
from config import * #(this line will give me access to all directory variables)
import matplotlib.pyplot as plt

def data_gen():

    KS_M = np.zeros((3,9))

    for k in range(0, 2+1):
        for j in range(0, 8+1):



        clusters = ['a0370','a1689','a1835','a2219','cl0024', 'ms0451','ms1054','rxj0152','rxj1347']
        clusname = cluster[j]

        bands = ['PSW', 'PMW', 'PLW']
        band = bands[k]

        I_0 = np.zeros((h+1))
        xbin = np.zeros((l+1))






def plot_hist() :
    b = 0
    dI = abs((max(I_0)-min(I_0))/l) # this is the binsize
    for b in range(len(l)-1) :
        xbin[b] = min(I_0)+ (b+0.5)*dI
        b += 1

    hist, binedges = np.histogram(d '''Some kind of data''')
                    # Default uses 10 equally sized bins and returns
                    # a tuple, includes one more bin edge than there are bins
                    #def histogram(a, bins=10, range=None, normed=False, weights=None,density=None):

    data_file = '/sz/szout_' + clusname + '.sav'
    data = scipy.io.readsav(CLUSDATA + data_file,python_dict = True)



    n, bins, patches = plt.hist(x=d, bins='auto',
                                alpha=0.7, rwidth=0.85)


    plt.xlabel('I_0 (mJy/sr)')
            #may need to fix notation for this text
    plt.ylabel('Probability')
    plt.title(string('I_0 for' + clusname)')
                    #should make the title based on the cluster name

    #Below is from the article but not clear what it is doing
    # Set a clean upper y-axis limit.
    maxfreq = n.max()
    plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)



    # A clearer version from the matplotlib documentation
    n_bins = 20

    fig, axs = plt.subplots(1, 2)

    # We can set the number of bins with the `bins` kwarg
    axs[0].hist(x, bins=n_bins)
    axs[1].hist(y, bins=n_bins)
