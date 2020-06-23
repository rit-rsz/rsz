################################################################################
# NAME : clus_likelihood
# DATE STARTED : Febuary 25th, 2020
# AUTHORS :  Benjamin Vaughan, Victoria Butler
# PURPOSE : This code is meant to do the likelihood analysis for the y, Te values
# for a cluster given the dI of the SZe effect in Mjy / Sr with some uncertainty.
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import sys
sys.path.append('../utilities')
sys.path.append('source_handling')
sys.path.append('reduc')
sys.path.append('sz')
# from clus_sim_hist import *
import config
from scipy.stats import norm, chi2
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib
import numpy as np
from IB_model import *
from random import *
from PIL import Image
from clus_get_relsz import *
from clus_sz_grids import *
from scipy.stats import rv_discrete
from scipy.interpolate import griddata
from scipy.special import gamma
from scipy.optimize import curve_fit

def fit_func(x, a, b, c):
    return a * np.exp(-1 * (x - b)**2 / (2*c**2))

def delta_chi_square_test(data, data_sig, y_nom, t_nom, t_sig):
    bolo_dI, x_pos, spectrum, xaxis = run_szpack('BOLOCAM', y_nom, t_nom + t_sig)
    PMW_dI, x_pos, spectrum, xaxis = run_szpack('PMW', y_nom, t_nom + t_sig)
    PLW_dI, x_pos, spectrum, xaxis = run_szpack('PLW', y_nom, t_nom + t_sig)
    ref_bolo_dI, x_pos, spectrum, xaxis = run_szpack('BOLOCAM', y_nom, t_nom)
    ref_PMW_dI, x_pos, spectrum, xaxis = run_szpack('PMW', y_nom, t_nom)
    ref_PLW_dI, x_pos, spectrum, xaxis = run_szpack('PLW', y_nom, t_nom)

    model = [bolo_dI, PMW_dI, PLW_di]
    chi1 = chi_square_test(data,model,data_sig)

    model = [ref_bolo_dI, ref_PMW_dI, ref_PLW_dI]
    chi2 = chi_square_test(data,model,data_sig)

    diff = chi1 - chi2
    print(diff)



def thermal(v, I, y):
    '''
    This is for use in calculating the cannonical effect for over plotting on top of the output graph
    Inputs: v (float array) - dimensionless frequency (see Carlstrom 2012)
            I (float) - Intensity
            y (float) - the compton y parameter ( this should be y_0 not y_{R2500})
    Outputs: outputs a spectrum of the SZe effect at the given frequencies.
    '''
    func1 = (((v**4)*np.exp(v))/(np.exp(v)-1)**2) * ((v*((np.exp(v) + 1)/(np.exp(v) - 1)))-4)
    return func1*(I*y)

def confidence_levels(data):
    '''
    Here we calculate the 1 sigma confidence intervals of the histogram plots
    for y, and Te.
    Inputs: data (float array) - this is meant to be a 1D array of data points representing the likelihood of the parameter of interest
    Outputs: lower (int) - an integer value corresponding to the index of the lower confidence level
             upper (int) - an integer value corresponding to the index of the upper confidence level
    '''
    sum = 0
    arr = data.copy() #copy the array so we don't have pass by reference issues and set our data array to zero.
    while sum < 0.683: #0.683 is the decimal percent for 1 sigma.
        #find the maximum likelihood, add that value to a summer and then set the array at that index to zero.
        #repeat this process until the sum is at the confidence level that we want.
        max_x = np.where(arr == np.amax(arr))
        max = arr[max_x[0]]
        sum += max
        arr[max_x[0]] = 0
    #this array represents the locations of all the zeros where we found a maximum in the likelihood array.
    #assuming that all points are initially non-zero and that the graph is continous the edges of this array should correspond
    #to the confidence intervals
    conf_integral = np.where(arr == 0)[0]
    lower = conf_integral[0] #pick out the lower edge
    upper = conf_integral[-1] #pick out the upper edge
    return lower, upper

def sum_confidence(data):
    '''
    In this routine we want to calculate the confidence contours. We do that by iterativley looking through the array for
    the next highest maximum likelihood and adding it to a sum. Once that sum reaches a critical value (i.e 1 sigma, 2 sigma, or 3 sigma)
    we move on to the next confidence interval. In doing this process we set a 2D array of (M, M) to values of 1, 2, or 3 corresponding to
    the confidence level at that point in the 2D likelihood plot.
    Inputs: data (float array) - an (M,M) array of data points corresponding to a 2D chart of the likelihood for y, Te values.
    Outputs:  conf_arr (float array) - an (M, M) array of values at 1, 2, 3 corresponding to 1 sigma, 2sigma, or 3 sigma on the
    likelihood graph.
    '''
    sum = 0
    conf_arr = np.zeros(data.shape)
    while sum < 0.683: #1 sigma confidence level
        max_x, max_y = np.where(data == np.amax(data))
        max = data[max_x[0], max_y[0]]
        conf_arr[max_x[0], max_y[0]] = 3
        data[max_x[0], max_y[0]] = 0
        sum += max
    while sum < .955: #2sigma confidence level
        max_x, max_y = np.where(data == np.amax(data))
        max = data[max_x[0], max_y[0]]
        conf_arr[max_x[0], max_y[0]] = 2
        data[max_x[0], max_y[0]] = 0
        sum += max
    # while sum < .997: #3 sigma confidence level
    #     max_x, max_y = np.where(data == np.amax(data))
    #     max = data[max_x[0], max_y[0]]
    #     conf_arr[max_x[0], max_y[0]] = 1
    #     data[max_x[0], max_y[0]] = 0
    #     sum += max
    return conf_arr

def chi_square_test(data,model,sigma):
    #this is a standard chi_squared calculation
    chi_ele = []
    for i in range(len(data)):
        chi_ele.append((data[i] - model[i])**2 / sigma[i]**2)
    return np.sum(chi_ele)

def create_test_data():
    '''
    this is for testing, eventually we will integrate this into more of a pipeline and not have to pass in the arguments like this.
    all this script does is take the values for dI that we have computed thus far and put them into a format that the likelihood analysis
    script can use.
    '''

    bolo_sz = -0.940
    psw_sz = 0.11
    pmw_sz = 0.14
    plw_sz = 0.55
    test_data = [bolo_sz, pmw_sz, plw_sz]
    bolo_sig = 0.053
    psw_sig = 0.05
    pmw_sig = 0.04
    plw_sig = 0.08
    test_sigma = [bolo_sig, pmw_sig, plw_sig]
    return test_data, test_sigma

def clus_likelihood(data, sigma, name='rxj1347',samples=10):
    '''
    Inputs: data (list) - a list of the dI values for each band must be in order ('BOLO', 'PMW', 'PLW')
            sigma (list) - a list of the associated uncertainties for the above dI values
            name (str) - the name of the cluster of interest
            samples (int) (depreciated) - the length / width of the look up grid.
    Outputs: y, Te values with uncertainties and creates a corner contour plot of the uncertainties
             and a best fit plot for the most probable y, Te values with the cannonical spectrum at Te = 0 for reference
    '''


    if not os.path.isfile(config.OUTPUT + 'sz_grids/%s_sz_grid_BOLOCAM.npy' %(name)):
        print('No grids detected for grid search, please run clus_sz_grid.py')
    sz_grid_0 = np.load(config.OUTPUT + 'sz_grids/%s_sz_grid_BOLO.npy' % (name), allow_pickle=True)
    sz_grid_1 = np.load(config.OUTPUT + 'sz_grids/%s_sz_grid_PSW.npy' % (name), allow_pickle=True)
    sz_grid_2 = np.load(config.OUTPUT + 'sz_grids/%s_sz_grid_PMW.npy' % (name), allow_pickle=True)
    sz_grid_3 = np.load(config.OUTPUT + 'sz_grids/%s_sz_grid_PLW.npy' % (name), allow_pickle=True)
    ys = np.load(config.OUTPUT + 'sz_grids/%s_ys_grid.npy' % name)
    ts = np.load(config.OUTPUT + 'sz_grids/%s_ts_grid.npy' % name)
    spec = np.load((config.OUTPUT + 'sz_grids/%s_spectrum_y.npy' % name), allow_pickle=True)
    xaxis = np.load((config.OUTPUT + 'sz_grids/%s_spectrum_x.npy' % name), allow_pickle=True)
    param_grid = np.load(config.OUTPUT + 'sz_grids/%s_y_t_grid.npy' % name, allow_pickle=True)
    #this is a place holder right now, ys and ts values are getting saved over the weekend.
    # ys = np.linspace(0, 12e-4, 5)
    # ts = np.linspace(0, 30, 5)

    #this is for calculating the bias once we have the simulations working currently not being used and not validated in anyway.

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
    # dI, x_pos, spectrum, xaxis = run_szpack('PLW', 0, ys[40], 0.3)

    #create data structures to hold the chi_square values and the likelihood values.
    like_li = np.zeros((len(ys),len(ts)))
    chi_grid = np.zeros((len(ys), len(ts)))

    df = len(data) - 1 #degrees of freedom
    for i in range(len(ys)):
        for j in range(len(ts)):
            #calculate the chisquare value for each y, Te combination
            chi_stat = chi_square_test(data, [sz_grid_0[i,j], sz_grid_2[i,j], sz_grid_3[i,j]], sigma)
            chi_grid[i,j] = chi_stat
            #calculate the likelihood from the chisquare value
            like_li[i,j] = ( chi_stat **(df/2 - 1) * np.exp(-1 * chi_stat/2) ) / (2 ** (df/2) * gamma(df / 2))

    #normalize the likelihood so that it sums to 1.
    like_norm = like_li / np.sum(like_li)

    y_max, t_max = np.where(like_norm == np.max(like_norm))
    print(param_grid[y_max, t_max, 0] * 0.163, param_grid[y_max, t_max, 1])
    #create data structures to hold y, Te likelihood values
    y_likelihood = np.zeros(len(ys))
    t_likelihood = np.zeros(len(ts))

    #calculate the marginalized likelihood for the y, Te values.
    for i in range(len(ys)):
        y_likelihood = np.add(y_likelihood, like_norm[:, i])
        t_likelihood = np.add(t_likelihood, like_norm[i, :])

    #pick out the most probable y and Te values
    y_max = np.where(y_likelihood == np.max(y_likelihood))
    t_max = np.where(t_likelihood == np.max(t_likelihood))


    #find the confidence levels.
    low_y, up_y = confidence_levels(y_likelihood)
    low_t, up_t = confidence_levels(t_likelihood)


    #convert to y_R2500
    ys = np.multiply(ys, 0.163)

    #create a meshgrid for contour / color map plots
    X, Y = np.meshgrid(ys, ts)


    matplotlib.rcParams.update({'font.size': 8})

    #begin the process of creating plots !
    fig, axs = plt.subplots(2, 2)

    #y = 1.44e-4 (-0.065,+0.073), Te = 8.3 (-2.7,+2.9)
    #print the best fits to screen !
    print('best fit y: %.3E - %.3E / + %.3E' % (ys[y_max], ys[y_max] - ys[low_y], ys[up_y] - ys[y_max]))
    print('best fit t: %.3E - %.3E / + %.3E' % (ts[t_max], ts[t_max] - ts[low_t], ts[up_t] - ts[t_max]))

    # p, e = curve_fit(fit_func, ys, y_likelihood, p0=[0.0140, 1.44e-4, 1e-5])
    # line = fit_func(ys, *p)
    # delta_chi_square_test(data, sigma, ys[y_max], ts[t_max], ts[up_t])


    #plot the y_2500 histogram
    # axs[0,0].plot(ys, line)
    axs[0,0].plot(ys, y_likelihood)
    axs[0,0].axvline(ys[y_max], linestyle='dashed')
    axs[0,0].axvline(ys[low_y], linestyle='dashed')
    axs[0,0].axvline(ys[up_y], linestyle='dashed')
    axs[0,0].set_yticklabels([])
    axs[0,0].set_xticklabels([])
    # plt.savefig(config.OUTPUT + 'y_T_fits/%s_ys_hist.png'%name)
    # plt.clf()

    # p, e = curve_fit(fit_func, ts, t_likelihood, p0=[0.0140, 8.5, 3])
    # line = fit_func(ts, *p)
    #plot the Te histogram
    # axs[1,1].plot(line, ts)
    axs[1,1].plot(t_likelihood, ts)
    axs[1,1].axhline(ts[t_max], linestyle='dashed')
    axs[1,1].axhline(ts[low_t], linestyle='dashed')
    axs[1,1].axhline(ts[up_t], linestyle='dashed')
    axs[1,1].set_yticklabels([])
    axs[1,1].set_xticklabels([])

    # axs[1,0].savefig(config.OUTPUT + 'y_T_fits/%s_ts_hist.png'%name)
    # plt.clf()

    #we have to switch axes for plotting the likelihood maps
    like_norm = np.swapaxes(like_norm, 0, 1)
    chi_grid = np.swapaxes(chi_grid, 0, 1)

    #plot the contour plot
    axs[1,0].pcolormesh(X, Y, like_norm) #for displaying the underlying pdf
    axs[1,0].set_xlabel('$<y>_{R2500}$')
    axs[1,0].set_ylabel('$T_e$ [KeV]')
    axs[1,0].ticklabel_format(axis='x', style='sci', scilimits=(-4,-4))

    # axs[1,0].colorbar().set_label('Normalized Likelihood')
    conf_interval = sum_confidence(like_norm)
    axs[1,0].contour(ys, ts, conf_interval, colors='gray')
    axs[1,0].scatter(ys[y_max], ts[t_max], marker='x', color='black')
    axs[1,0].scatter(1.29e-4, 9.47, marker='^', color='black')


    #fixing layout
    fig.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.05, hspace=0.05)
    # plt.tight_layout()
    fig.delaxes(axs[0,1])
    fig.savefig(config.OUTPUT + 'y_T_fits/%s_y_T_analysis.png' % name)
    plt.clf()


    #below is more math to create the cannonical SZ spec at Te = 0
    T_0     = 2.725                   #CMB temperature, K
    k_B     = 1.3806503e-23           #Boltzmann constant, J/K
    h       = 6.626068e-34            #Planck constant, J s
    c       = 2.99792458e8            #m/s
    GHztoHz = 1.0e9                   #GHz -> Hz
    I = 2*((k_B*T_0)**3/(h*c)**2) * 10**(26) * 10**(-6)

    v = (h/(k_B * T_0))*np.linspace(1*10**(9),1700*10**9,100)

    #create the cannonical spectrum .
    cannon_spec = thermal(v, I, ys[y_max] / 0.163)

    #convert from dimnesionless frequency to GHz.
    cannon_x = np.multiply( T_0 * k_B / (h * 1e9), v)
    best_x = np.multiply(T_0 * k_B / (h * 1e9), xaxis[y_max, t_max].flatten())

    #plot the best fitting spectrum with data points and the cannonical spectrum.
    # plt.plot(best_x,spec[y_max,t_max].flatten(), label='best_fit')
    # x_vals = [140.9, 856.55, 599.585]
    # plt.xlim(0, 1000)
    # plt.errorbar(x_vals, data, yerr=sigma, label='data', linestyle='None', marker='o',markersize=4, ecolor='black', markerfacecolor='black', markeredgecolor='black')
    # plt.plot(cannon_x, cannon_spec, linestyle='dashed', label='cannonical')
    # plt.legend(loc='best')
    # plt.xlabel('$\\nu$ (GHz)')
    # plt.ylabel('$\Delta$ $I_0$ [MJy/sr]')
    # # plt.set_position([0.7, 0.7, 0.95,0.95])
    # plt.savefig(config.OUTPUT + 'y_T_fits/%s_max_spectrum.png' % name)
    # plt.clf()
    # make_comp_image(name)




if __name__ == '__main__' :

    data, sigma = create_test_data()
    clus_likelihood(data, sigma, name='rxj1347',samples=100)
