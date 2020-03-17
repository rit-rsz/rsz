################################################################################
# NAME : test.py
# DATE STARTED : October 25, 2019
# AUTHORS : Benjamin Vaughan
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

def chi_square_test(data,model,sigma):
    total_chi = []
    for i in range(3):
        total_chi.append((data[i] - model[i])**2 / (2*sigma[i]**2))
    return np.sum(total_chi)

def plot_cov_ellipse(points, center, conf=0.683, color='orange'):

    # q = 2 * norm.cdf(1) - 1
    # r2 = chi2.ppf(q, 2)
    cov = np.cov(points, rowvar=False)
    val, vec = np.linalg.eigh(cov)
    width = 2 * sqrt(val[:,None][0]*conf)
    height = 2 * sqrt(val[:,None][1]*conf)
    rotation = np.degrees(np.arctan2(*vec[::-1, 0]))
    # Width and height are "full" widths, not radius
    ellip = Ellipse(xy=center, width=width, height=height, angle=rotation, facecolor='none',edgecolor=color, lw = 3, linestyle='--')

    return ellip

def clus_likelihood(nsim=0,name='rxj1347',samples=100,step=0.001):

    #yt_grid and sz_grid should be indexed the same
    if os.path.isfile(config.HOME + 'outputs/%s_sz_grid.npy' %(name)):
        sz_grid = np.load(config.HOME + 'outputs/%s_sz_grid.npy' %(name),allow_pickle=True)
        input_yt = np.load(config.HOME + 'outputs/%s_input_grid.npy' %(name),allow_pickle=True)
        param_grid = np.load(config.HOME + 'outputs/%s_param_grid.npy' %(name))
        ys = np.load(config.HOME + 'outputs/%s_y.npy' %(name))
        ts = np.load(config.HOME + 'outputs/%s_t.npy' %(name))
    else :
        sz_grid, input_yt, param_grid, ys, ts = clus_szgf(nsim,name,samples,step) # input_yt should have one DI for each band

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

    # what is the critical p value ?
    crit_p = chi2.ppf(q = 0.683, df = 2)
    print('critical value 68.3% confidence :',crit_p)

    # grab real sz fit params
    # sz_fit_real = np.load(config.HOME + 'outputs/%s_fit.npy'%(name))
    # print(sz_fit_real)

    # retreive and separate all dI for each band and y/t pair
    # subtract bias for each band
    ''' TESTING until we get realistic SZ values and sigma'''
    like_li = []
    sz_grid_0 = [x.get(0) for x in sz_grid]
    sz_grid_1 = [x.get(1) for x in sz_grid]
    sz_grid_2 = [x.get(2) for x in sz_grid]
    
    print(sz_grid_2[0:10])

    x = np.arange(0, 100)
    y = np.arange(0, 100)
    y = y[:, np.newaxis]
    rc = 15.0
    beta = 0.8
    beta_map = (1. + (np.sqrt((x-50)**2 + (y-50.0)**2) / rc)**2)**((1. - 3. * beta) / 2.)
    plt.imshow(beta_map)
    plt.colorbar()
    plt.savefig('beta_map.png')
    ''' #################################################################################################################### '''
    # make chi square test for each final sz amplitude
    # use sz_fit_real for real DIs
    seed(102793)
    for i in range(len(sz_grid)):
        # chi_stat = chi_square_test([sz_fit_real[0][0],sz_fit_real[1][0],sz_fit_real[2][0]],
        #                             [sz_grid[i].get(0),sz_grid[i].get(1),sz_grid[i].get(2)],
        #                             [bias[0],bias[1],bias[2]])
        chi_stat = chi_square_test([sz_grid[i].get(0)-uniform(0,0.001),sz_grid[i].get(1)-uniform(0,0.001),sz_grid[i].get(2)-uniform(0,0.001)],
                                    [sz_grid[i].get(0),sz_grid[i].get(1),sz_grid[i].get(2)],
                                    [0.05,0.05,0.05])

        like_li.append(np.exp(chi_stat/-2.0))

    like_norm = list(like_li / np.sum(like_li))
    # like_norm = list(like_norm * beta_map.flatten())
    like_max = like_norm.index(max(like_norm))

    # compute likelihood function
    DI = np.array(like_norm).reshape(len(ys),len(ts))
    center = [param_grid[like_max][0],param_grid[like_max][1]]
    YS,TS = np.meshgrid(ys,ts)
    ellip = plot_cov_ellipse(param_grid,center)
    ellip2 = plot_cov_ellipse(param_grid,center,color='green',conf=0.85)
    ellip3 = plot_cov_ellipse(param_grid,center,color='red',conf=0.95)
    fig,ax = plt.subplots()
    plt.pcolormesh(YS,TS,DI)
    ax.add_artist(ellip)
    ax.add_artist(ellip2)
    ax.add_artist(ellip3)
    plt.colorbar().set_label('Likelihood')
    ells = [ellip,ellip2,ellip3]
    ax.legend(ells,['68.3%','85%','95%'])

    # compute the 68.3% confidence contours
    # cs = plt.contour(YS,TS,DI,levels=[0.683],colors=('k',),linestyles=('-',),linewidths=(2,))
    # plt.clabel(cs,fmt='%.3f',colors='k',fontsize=14)
    plt.xlabel('Compton Y')
    plt.ylabel('Temperature [K]')
    plt.savefig('likelihood.png')
    plt.clf()

if __name__ == '__main__' :
    clus_likelihood(nsim=0,name='rxj1347')
