################################################################################
# NAME : fitsz.py
# DATE STARTED : June 24, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE :
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import numpy as np
import matplotlib.pyplot as plt
import sys
from math import *
sys.path.append('../utilities')
from get_clusparams import *
from compute_rings import *
from get_data import *
from scipy.stats import linregress
from astropy.convolution import convolve_fft
from scipy.interpolate import splrep
from gaussian import *

def fitsz(radave, params, beam=None, maxlim=3600, minlim=0, noweight=1, superplot=1, verbose=1):

    # init params
    ncols = len(radave)
    increment = np.zeros(ncols)
    offsets = np.zeros(ncols)
    fit = np.zeros((2, ncols))

    for i in range(ncols):
        print(params['rc'])
        print(radave[i]['midbin'])
        print(radave[i]['midbin'].shape)
        rc = np.tile(params['rc'], radave[i]['midbin'].shape)
        beta = np.tile(params['beta'], radave[i]['midbin'].shape)
        if  not beam:
            roft = (1. + (radave[i]['midbin'] / rc)**2)**((1. - 3. * beta / 2.))
            plt.plot(roft)
            plt.show()
        #     for j in range(radave[i]['midbin'].shape[0]):
        #         for k in range(radave[i]['midbin'].shape[1]):
        #             roft = (1.0 + (radave[i]['midbin'] / rc)**2)**((1.0 - 3.0 * beta) / 2.0)
        else:
            print('this part is not written yet')
            #this stuff is commented out because I don't understand it also it's probably wrong.
            print(np.max(radave[i]['midbin']))
            length = int(np.max((radave[i]['midbin'])))
            xrad = np.arange(-1 * length, length)
            # xrady = np.arange(0, int(np.max(radave[i]['midbin'])))
            # xrady = xrady[:, np.newaxis]
            # xrad = np.array([xrad1, xrad2])

            # xrad = np.array([xrad1, xrad2])
            xrc = np.tile(params['rc'],2 * length)
            xbeta = np.tile(params['beta'],2 * length)
            roftprime = (1. + (xrad / xrc)**2)**((1. - 3. * xbeta) / 2.)
            pf = -1 * (2.18 * log(2) / beam[i]**2)
            psf = np.exp(pf * xrad**2)
            psf = psf / np.sum(psf)
            plt.plot(psf)
            plt.show()
            roftprimep = convolve_fft(roftprime, psf)
            plt.plot(roftprimep)
            plt.show()
            roftprimep = roftprimep / np.max(np.abs(roftprimep))
            plt.plot(roftprimep)
            plt.show()
            # roft = splrep(y=roftprimep, x=xrad, xe=radave[i]['midbin'])

        if superplot:
            plt.plot(radave[i]['midbin'], radave[i]['fluxbin'])
            plt.show()

        ind = []
        for j in range(len(roft)):
            if np.isfinite(radave[i]['fluxbin'][j]) and radave[i]['midbin'][j] <= maxlim:
                ind.append(j)
        ind = np.asarray(ind)

        slope, intercept, r_value, p_value, std_err = linregress(roft[ind], radave[i]['fluxbin'][ind])

        A = slope
        B = intercept

        print(A, B)

        if superplot:
            y = slope * roft + intercept
            plt.plot(roft, radave[i]['fluxbin'], c='blue')
            plt.plot(roft, y, c='red')
            plt.show()
            #plot the new stuffs, again could probably be done in matplotlib, but i need to do real testing to find out.

        increment[i] = slope
        offsets[i] = intercept
        fit[:, i] = np.array([slope, intercept])

    return fit, offsets, increment

if __name__ == '__main__':
    midbin = np.array([  15.09459459,   45.28378378,   75.47297297,  105.66216216,
        135.85135135,  166.04054054,  196.22972973,  226.41891892,
        256.60810811,  286.7972973 ,  316.98648649,  347.17567568,
        377.36486486,  407.55405405,  437.74324324,  467.93243243,
        498.12162162,  528.31081081,  558.5       ,  588.68918919,
        618.87837838,  649.06756757,  679.25675676,  709.44594595,
        739.63513514,  769.82432432,  800.01351351,  830.2027027 ,
        860.39189189,  890.58108108,  920.77027027,  950.95945946,
        981.14864865, 1011.33783784, 1041.52702703, 1071.71621622,
       1101.90540541])
    radbin = np.array([   0.        ,   30.18918919,   60.37837838,   90.56756757,
        120.75675676,  150.94594595,  181.13513514,  211.32432432,
        241.51351351,  271.7027027 ,  301.89189189,  332.08108108,
        362.27027027,  392.45945946,  422.64864865,  452.83783784,
        483.02702703,  513.21621622,  543.40540541,  573.59459459,
        603.78378378,  633.97297297,  664.16216216,  694.35135135,
        724.54054054,  754.72972973,  784.91891892,  815.10810811,
        845.2972973 ,  875.48648649,  905.67567568,  935.86486486,
        966.05405405,  996.24324324, 1026.43243243, 1056.62162162,
       1086.81081081, 1117.        ])
    radave = [[],[]]
    pmw = {'midbin' : midbin,
           'radbin' : radbin}
    midbin = np.array([  15.2972973 ,   45.89189189,   76.48648649,  107.08108108,
        137.67567568,  168.27027027,  198.86486486,  229.45945946,
        260.05405405,  290.64864865,  321.24324324,  351.83783784,
        382.43243243,  413.02702703,  443.62162162,  474.21621622,
        504.81081081,  535.40540541,  566.        ,  596.59459459,
        627.18918919,  657.78378378,  688.37837838,  718.97297297,
        749.56756757,  780.16216216,  810.75675676,  841.35135135,
        871.94594595,  902.54054054,  933.13513514,  963.72972973,
        994.32432432, 1024.91891892, 1055.51351351, 1086.10810811,
       1116.7027027 ])
    radbin = np.array([   0.        ,   30.59459459,   61.18918919,   91.78378378,
        122.37837838,  152.97297297,  183.56756757,  214.16216216,
        244.75675676,  275.35135135,  305.94594595,  336.54054054,
        367.13513514,  397.72972973,  428.32432432,  458.91891892,
        489.51351351,  520.10810811,  550.7027027 ,  581.2972973 ,
        611.89189189,  642.48648649,  673.08108108,  703.67567568,
        734.27027027,  764.86486486,  795.45945946,  826.05405405,
        856.64864865,  887.24324324,  917.83783784,  948.43243243,
        979.02702703, 1009.62162162, 1040.21621622, 1070.81081081,
       1101.40540541, 1132.        ])
    plw = {'midbin' : midbin,
           'radbin' : radbin}
    radave = [[],[]]
    radave[0] = pmw
    radave[1] = plw
    beam = [18, 25, 36]
    params, err = get_clus_params('a0370')
    maps, err = get_data('a0370')
    psw = fits.open('../fits_files/xid_9_subtracted_a0370_PSW.fits')
    pmw = fits.open('../fits_files/xid_9_subtracted_a0370_PMW.fits')
    plw = fits.open('../fits_files/xid_9_subtracted_a0370_PLW.fits')
    maps[0]['srcrm'] = psw[0].data
    maps[1]['srcrm'] = pmw[0].data
    maps[2]['srcrm'] = plw[0].data
    # radave = compute_rings(maps,params,30.0)
    print(radave[0]['midbin'])
    fit = fitsz(radave, params)
