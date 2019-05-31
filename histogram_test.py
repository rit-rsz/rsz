
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
import math

def data_gen(nsims,len):

    KS_M = np.zeros((3,9))
    cluster = ['a0370','a1689','a1835','a2219','cl0024', 'ms0451','ms1054','rxj0152','rxj1347']
    bands = ['PSW', 'PMW', 'PLW']

    for k in range(len(bands)-1):
        I_0_sim = np.zeros((nsims+1))
        xbin = np.zeros((len+1))
        band = bands[k]

        for j in range(len(clusters)-1):
            clusname = cluster[j]

            n = 0
            for i in range(200,200+nsims-1) :
                data_file = 'sz/szout_0' i + '_' + clusname + '.sav'
                data = scipy.io.readsav(CLUSDATA + data_file,python_dict = True)
                arr_data = data.values()[0][k][0] # only works on python 2, returns increment
                ''' for python 3 use arr_data = list(data.values())[0][k][0] '''
                # filter out nan values if they exist #
                # final_data = [value for value in arr_data if not math.isnan(value)]
                I_0_sim[n] = arr_data
                n += 1
                plot_hist(j,k,len)



def plot_hist(j,k,l) :
'''It does make len(len) confusing though'''
''' I wondered about that, probably we should yeah'''
    b = 0
    dI = abs((max(I_0)-min(I_0))/l) # this is the binsize
    for b in range(len(l)-1) :
        xbin[b] = min(I_0)+ (b+0.5)*dI
        b += 1

    data_file = 'sz/szout_' + clusname + '.sav'
    data = scipy.io.readsav(CLUSDATA + data_file,python_dict = True)
    ########################## Data Array Structure ###########################
    # [increment,offset,band,radbin,midbin,fluxbin,errbin,rc,beta,clusname]
    ###########################################################################
    inc = data.values()[0][k][0] # only works on python 2, returns increment
    ''' for python 3 use arr_data = list(data.values())[0][k][0] '''
    I_0 = inc
#     We need to determine where this then needs to go

    xmin = min([I_0,min(xbin)])
    xmax = max([I_0,max(xbin)])
    n, bins, patches = plt.hist(x=arr_data, bins='auto', rwidth=0.85)
                                # rwidth is the width of the bars as a fraction of the bin width
    Gmax = float(abs(max(n)))
    N = len(xbin)
    Area = np.sum(G)
    Area_av = Area/max(G)

    plt.xlabel('I_0 (mJy/sr)')
    plt.ylabel('Probability')
    plt.title(string('I_0 for' + clusname)')
    plt.show()


'''    hist, binedges = np.histogram("Some kind of data")
                    # Default uses 10 equally sized bins and returns
                    # a tuple, includes one more bin edge than there are bins
                    #def histogram(a, bins=10, range=None, normed=False, weights=None,density=None):

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
'''

if __name__ == '__main__':
    data_gen(100,10)
