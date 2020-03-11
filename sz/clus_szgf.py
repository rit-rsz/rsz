################################################################################
# NAME : clus_szgf.py
# DATE STARTED : January 23, 2020
# AUTHORS : Victoria Butler
# PURPOSE : This produces an SZ amplitude for a given grid of y and Te,
#           then performs a Chi^2 fit test to determine the optimal parameter space.
# INPUTS :
#        nu  - the range of frequncies requested [fltarr, GHz]
#        y   - the Compton y parameter, usually O(10^-4) [float, unitless]
#        T_e - the electron temperature of the gas [float, keV]
#        vpec = cluster peculiar velocity, usually 0.0 [float, m/s]
#
# OUTPUTS : rSZ intensity at given SPIRE band frequency [MJy/sr]
# NOTES :
#
# REVISION HISTORY :
################################################################################

import numpy as np
from math import *
import os, sys, time
sys.path.append('../utilities')
from clus_get_lambdas import *
from config import *
sys.path.append('../sz')
sys.path.append('../source_handling')
from clus_get_data import *
from clus_get_relsz import *
import multiprocessing as mp

def clus_szgf(nsim,name,samples,step):
    yin_coeff = [4.56,8.5,8.7,4.45,7.43,8.96,0.93,7.42,5.63,3.53,3.72,14.08]
    f_yin = [x*1e-4 for x in yin_coeff]

    f_tin = [7.2,10.1,7.65,6.7,9.81,9.16,4.5,8.62,7.8,7.2,5.5,10.88]

    clusters = ['a0370','a1689','a1835','a2218','a2219','a2390',
              'cl0024','ms0451','ms1054','ms1358','rxj0152','rxj1347']

    tin = f_tin[clusters.index(name)]
    yin = f_yin[clusters.index(name)]

    # for i in range(len(clusters)):
    maps,err = clus_get_data(name,nsim)
    ys = np.linspace(yin - (step*samples/2.0*1e-4), yin + (step*samples/2.0*1e-4), samples) # 100 samples, 0.01 step
    ts = np.linspace(tin - (step*samples/2.0), tin + (step*samples/2.0), samples)
    param_grid = np.array(np.meshgrid(ys,ts)).T.reshape(-1,2) #[0,0] = ys,ts
    a = time.time()
    master_sz = []

    for j in range(len(param_grid)):
        y_c = param_grid[j][0]
        t_c = param_grid[j][1]

        # default is to have all of the 3 maps returned
        manager = mp.Manager()
        sz_amp = manager.dict()
        p1 = mp.Process(target=run_szpack, args=(maps,sz_amp,0,nsim,y_c,t_c))
        p2 = mp.Process(target=run_szpack, args=(maps,sz_amp,1,nsim,y_c,t_c))
        p3 = mp.Process(target=run_szpack, args=(maps,sz_amp,2,nsim,y_c,t_c))
        p1.start()
        p2.start()
        p3.start()
        p1.join()
        p2.join()
        p3.join()
        master_sz.append(sz_amp.copy()) # have to return copy since sz_amp is DictProxy Object

    # this is to calculate the central input value used for average
    master_yt = []
    mngr = mp.Manager()
    input_yt = mngr.dict()
    p1 = mp.Process(target=run_szpack, args=(maps,input_yt,0,nsim,yin,tin))
    p2 = mp.Process(target=run_szpack, args=(maps,input_yt,1,nsim,yin,tin))
    p3 = mp.Process(target=run_szpack, args=(maps,input_yt,2,nsim,yin,tin))
    p1.start()
    p2.start()
    p3.start()
    p1.join()
    p2.join()
    p3.join()
    master_yt.append(input_yt.copy())

    b = time.time()
    print('TIME ELAPSED:' , b-a)
    np.save(config.HOME + '/outputs/%s_sz_grid.npy'%(name),master_sz,allow_pickle=True)
    np.save(config.HOME + '/outputs/%s_input_grid.npy'%(name),master_yt,allow_pickle=True)
    np.save(config.HOME + '/outputs/%s_y.npy'%(name),ys)
    np.save(config.HOME + '/outputs/%s_t.npy'%(name),ts)

    return master_sz, master_yt, ys, ts

def run_szpack(maps,sz_amp,band,isim,y_c,t_c):
    nu = 3e5 / clus_get_lambdas((maps[band]['band']))
    dI,errmsg = clus_get_relsz(isim,nu,band,y=y_c,te=t_c) # dI = [MJy/sr]
    szin = dI / maps[band]['calfac'] # converts to Jy/beam
    sz_amp[band] = szin
    return

if __name__ == '__main__':
    clus_szgf()
