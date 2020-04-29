################################################################################
# NAME : test.py
# DATE STARTED : October 25, 2019
# AUTHORS : Victoria Butler , Benjamin Vaughan
# PURPOSE : This code is just for a first test run of actually running through
# the whole data analysis pipeline and returning graphs of the results.
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import sys
sys.path.append('utilities')
sys.path.append('source_handling')
sys.path.append('reduc')
sys.path.append('sz')
from clus_szgf import *
from clus_sim_hist import *
import config
from scipy.stats import norm, chi2
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
from IB_model import *
from random import *

def sum_confidence(data):
    sum = 0
    conf_arr = np.zeros(data.shape)
    while sum < .683:
        max_x, max_y = np.where(arr == np.amax(data))
        max = data[max_x[0], max_y[0]]
        conf_arr[max_x[0], max_y[0]] = 3
        data[max_x[0], max_y[0]] = 0
        sum += max
    while sum < .955:
        max_x, max_y = np.where(data == np.amax(data))
        max = data[max_x[0], max_y[0]]
        conf_arr[max_x[0], max_y[0]] = 2
        data[max_x[0], max_y[0]] = 0
        sum += max
    while sum < .997:
        max_x, max_y = np.where(data == np.amax(data))
        max = data[max_x[0], max_y[0]]
        conf_arr[max_x[0], max_y[0]] = 1
        data[max_x[0], max_y[0]] = 0
        sum += max
    return conf_arr

def chi_square_test(data,model,sigma):
    total_chi = []
    for i in range(3):
        total_chi.append((data[i] - model[i])**2 / (sigma[i]**2))
    return np.sum(total_chi)/2.0

def clus_likelihood(nsim=0,name='rxj1347',samples=10):

    # yt_grid and sz_grid should be indexed the same
    if os.path.isfile(config.HOME + 'outputs/%s_sz_grid.npy' %(name)):
        sz_grid_0 = np.load(config.HOME + 'outputs/%s_sz_grid_0.npy' %(name),allow_pickle=True)
        sz_grid_1 = np.load(config.HOME + 'outputs/%s_sz_grid_1.npy' %(name),allow_pickle=True)
        sz_grid_2 = np.load(config.HOME + 'outputs/%s_sz_grid_2.npy' %(name),allow_pickle=True)
        input_yt = np.load(config.HOME + 'outputs/%s_input_grid.npy' %(name),allow_pickle=True)
        param_grid = np.load(config.HOME + 'outputs/%s_param_grid.npy' %(name))
        ys = np.load(config.HOME + 'outputs/%s_y.npy' %(name))
        ts = np.load(config.HOME + 'outputs/%s_t.npy' %(name))
        master_x_0 = np.load(config.HOME + 'outputs/%s_master_x_0.npy'%(name),allow_pickle=True)
        master_x_1 = np.load(config.HOME + 'outputs/%s_master_x_1.npy'%(name),allow_pickle=True)
        master_x_2 = np.load(config.HOME + 'outputs/%s_master_x_2.npy'%(name),allow_pickle=True)
        master_y = np.load(config.HOME + 'outputs/%s_master_y.npy'%(name),allow_pickle=True)
        master_array = np.load(config.HOME + 'outputs/%s_x_array.npy'%(name),allow_pickle=True)
        terp_x_0 = np.load(config.HOME + 'outputs/%s_spectrum_x_0.npy'%(name),allow_pickle=True)
        terp_x_1 = np.load(config.HOME + 'outputs/%s_spectrum_x_1.npy'%(name),allow_pickle=True)
        terp_x_2 = np.load(config.HOME + 'outputs/%s_spectrum_x_2.npy'%(name),allow_pickle=True)
        terp_y = np.load(config.HOME + 'outputs/%s_spectrum_y.npy'%(name),allow_pickle=True)
        xaxis = np.load(config.HOME + 'outputs/%s_x_array.npy'%(name),allow_pickle=True)
    else :
        sz_grid_0, sz_grid_1, sz_grid_2, terp_x_0, terp_x_1, terp_x_2, input_yt, ys, ts, real_terp_x, real_terp_y, xaxis, real_x_0, real_x_1, real_x_2, terp_y = clus_szgf(nsim,name,samples) # input_yt should have one DI for each band

    # avg_dI = clus_sim_hist(nsim,name)
    # print('avg dI : ',avg_dI)
    # print('input dI :',input_yt[0])
    #
    # # calculate sz bias from pipeline
    # bias = [0]*3
    # for i in range(3):
    #     if avg_dI[i] < 0 :
    #         bias[i] = input_yt[0].get(i) + abs(avg_dI[i])
    #     else :
    #         bias[i] = avg_dI[i] - input_yt[0].get(i)
    # print('bias in dI : ',bias[0],bias[1],bias[2])

    # save the sz spectrum output for this particular realization
    for i in range(len(sz_grid_0)):
        for j in range(len(sz_grid_0)):
            plt.plot(xaxis[i,j],terp_y[i,j])
            plt.scatter(terp_x_0[i,j],sz_grid_0[i,j],color='orange',label='500 Micron')
            plt.scatter(terp_x_1[i,j],sz_grid_1[i,j],color='green',label='350 Micron')
            plt.scatter(terp_x_2[i,j],sz_grid_2[i,j],color='purple',label='250 Micron')
            plt.title('4 Band dI Fit for y: %.6f t: %.6f' %(ys[i],ts[j]))
            plt.xlabel('Dimensionless Frequency')
            plt.ylabel('$\Delta$I_0 [MJy/sr]')
            plt.legend()
            plt.ylim((-1.0,1.25))
            plt.savefig('sz_spectrum_%s_%s.png' %(i,j))
            plt.clf()

    # what is the critical p value ?
    crit_p = chi2.ppf(q = 0.683, df = 2)
    print('critical value 68.3% confidence :',crit_p)

    # grab real sz fit params
    # sz_fit_real = np.load(config.HOME + 'outputs/%s_fit.npy'%(name))
    # print(sz_fit_real)

    # retreive and separate all dI for each band and y/t pair
    # subtract bias for each band

    # make chi square test for each final sz amplitude
    # use sz_fit_real for real DIs

    # bias = 3*[0]
    # PSW : FWHM = 18.0 "/pix
    # Jy/beam -> MJy/sr = 115.888
    # MJy/sr -> Jy/beam = 86.29E-4
    # PMW : FWHM = 25.0 "/pix
    # Jy/beam -> MJy/sr = 60.076
    # MJy/sr -> Jy/beam = 16.65E-3
    # PLW : FWHM = 36.0 "/pix
    # Jy/beam -> MJy/sr = 28.972
    # MJy/sr -> Jy/beam = 34.52E-3
    maps, err = clus_get_data(name,nsim,verbose=1,sgen=None,nsim=nsim, testflag=0)

    # bias[0] = 0.0058 * maps[0]['calfac']
    # bias[1] = 0.0063 * maps[1]['calfac']
    # bias[2] = 0.0068 * maps[2]['calfac']
    # bias = [0.0058,0.0063,0.0068]

    cen = len(sz_grid_0) - 1

    like_li = np.zeros((samples,samples))
    for i in range(samples) :
        for j in range(samples) :
        # chi_stat = chi_square_test([sz_fit_real[0][0],sz_fit_real[1][0],sz_fit_real[2][0]],
        #                             [sz_grid[i].get(0),sz_grid[i].get(1),sz_grid[i].get(2)],
        #                             [bias[0],bias[1],bias[2]])
            chi_stat = chi_square_test([sz_grid_0[cen,cen],sz_grid_1[cen,cen],sz_grid_2[cen,cen]],
                                        [sz_grid_0[i,j],sz_grid_1[i,j],sz_grid_2[i,j]],
                                        [0.05,0.05,0.05])

            like_li[i,j] = np.exp(chi_stat/-2.0)

    like_norm = like_li / np.sum(like_li)
    max_x, max_y = np.where(like_norm == np.amax(like_norm))
    # like_max = like_norm[max_x[0],max_y[0]]
    # compute likelihood function

    center = [ys[max_x][0],ts[max_y][0]]
    fig,ax = plt.subplots()
    plt.pcolormesh(ys,ts,like_norm)
    ax.scatter(center[0],center[1])
    # conf = sum_confidence()
    plt.colorbar().set_label('Normalized Likelihood')

    # compute the 68.3% confidence contours
    plt.xlabel('Compton Y')
    plt.ylabel('Temperature [K]')
    plt.savefig('likelihood.png')
    plt.clf()

if __name__ == '__main__' :
    clus_likelihood(nsim=0,name='rxj1347')
