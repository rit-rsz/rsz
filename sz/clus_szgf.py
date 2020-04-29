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

def clus_szgf(nsim, name, samples):
    yin_coeff = [
     4.56, 8.5, 8.7, 4.45, 7.43, 8.96, 0.93, 7.42, 5.63, 3.53, 3.72, 4.08]
    f_yin = [x * 0.0001 for x in yin_coeff]
    f_tin = [7.2, 10.1, 7.65, 6.7, 9.81, 9.16, 4.5, 8.62, 7.8, 7.2, 5.5, 10.88]
    clusters = [
     'a0370', 'a1689', 'a1835', 'a2218', 'a2219', 'a2390',
     'cl0024', 'ms0451', 'ms1054', 'ms1358', 'rxj0152', 'rxj1347']
    tin = f_tin[clusters.index(name)]
    yin = f_yin[clusters.index(name)]
    maps, err = clus_get_data(name, nsim)
    ys = np.linspace(1e-05, (2.0 * yin), num=samples)
    ts = np.linspace(2.0, (2.0 * tin), num=samples)
    master_sz_0 = np.zeros((len(ys), len(ts)))
    master_sz_1 = np.zeros((len(ys), len(ts)))
    master_sz_2 = np.zeros((len(ys), len(ts)))
    spectrum_x_0 = np.zeros((len(ys), len(ts)))
    spectrum_x_1 = np.zeros((len(ys), len(ts)))
    spectrum_x_2 = np.zeros((len(ys), len(ts)))
    spectrum_y = np.zeros((len(ys), len(ts), samples ** 2))
    x_array = np.zeros((len(ys), len(ts), samples ** 2))
    a = time.time()
    for j in range(len(ys)):
        y_c = round(ys[j], 6)
        for k in range(len(ts)):
            t_c = round(ts[k], 6)
            manager = mp.Manager()
            sz_amp = manager.dict()
            terp_x = manager.dict()
            terp_y = manager.dict()
            xaxis = manager.dict()
            p1 = mp.Process(target=run_szpack, args=(maps, sz_amp, 0, nsim, y_c, t_c, terp_x, terp_y, xaxis))
            p2 = mp.Process(target=run_szpack, args=(maps, sz_amp, 1, nsim, y_c, t_c, terp_x, terp_y, xaxis))
            p3 = mp.Process(target=run_szpack, args=(maps, sz_amp, 2, nsim, y_c, t_c, terp_x, terp_y, xaxis))
            p1.start()
            p2.start()
            p3.start()
            p1.join()
            p2.join()
            p3.join()
            master_sz_0[(j, k)] = sz_amp.copy().get(0)
            master_sz_1[(j, k)] = sz_amp.copy().get(1)
            master_sz_2[(j, k)] = sz_amp.copy().get(2)
            spectrum_x_0[(j, k)] = terp_x.copy().get(0)
            spectrum_x_1[(j, k)] = terp_x.copy().get(1)
            spectrum_x_2[(j, k)] = terp_x.copy().get(2)
            spectrum_y[j, k, :] = terp_y.copy().get(0)
            x_array[j, k, :] = xaxis.copy().get(0)

    master_yt = []
    master_x_0 = []
    master_x_1 = []
    master_x_2 = []
    master_y = []
    master_array = []
    realx_array = []
    mngr = mp.Manager()
    input_yt = mngr.dict()
    terp_x = manager.dict()
    terp_y = manager.dict()
    xaxis = manager.dict()
    p1 = mp.Process(target=run_szpack, args=(maps, input_yt, 0, nsim, yin, tin, terp_x, terp_y, xaxis))
    p2 = mp.Process(target=run_szpack, args=(maps, input_yt, 1, nsim, yin, tin, terp_x, terp_y, xaxis))
    p3 = mp.Process(target=run_szpack, args=(maps, input_yt, 2, nsim, yin, tin, terp_x, terp_y, xaxis))
    p1.start()
    p2.start()
    p3.start()
    p1.join()
    p2.join()
    p3.join()
    master_yt.append(input_yt.copy())
    master_x_0.append(terp_x.copy().get(0))
    master_x_1.append(terp_x.copy().get(1))
    master_x_2.append(terp_x.copy().get(2))
    master_y.append(terp_y.copy().get(0))
    realx_array.append(xaxis.copy().get(0))
    b = time.time()
    print('TIME ELAPSED:', b - a)
    np.save((config.HOME + 'outputs/%s_sz_grid_0.npy' % name), master_sz_0, allow_pickle=True)
    np.save((config.HOME + 'outputs/%s_sz_grid_1.npy' % name), master_sz_1, allow_pickle=True)
    np.save((config.HOME + 'outputs/%s_sz_grid_2.npy' % name), master_sz_2, allow_pickle=True)
    np.save((config.HOME + 'outputs/%s_input_grid.npy' % name), master_yt, allow_pickle=True)
    np.save((config.HOME + 'outputs/%s_master_x_0.npy' % name), master_x_0, allow_pickle=True)
    np.save((config.HOME + 'outputs/%s_master_x_1.npy' % name), master_x_1, allow_pickle=True)
    np.save((config.HOME + 'outputs/%s_master_x_2.npy' % name), master_x_2, allow_pickle=True)
    np.save((config.HOME + 'outputs/%s_master_y.npy' % name), master_y, allow_pickle=True)
    np.save((config.HOME + 'outputs/%s_x_array.npy' % name), master_array, allow_pickle=True)
    np.save(config.HOME + 'outputs/%s_y.npy' % name, ys)
    np.save(config.HOME + 'outputs/%s_t.npy' % name, ts)
    np.save((config.HOME + 'outputs/%s_spectrum_x_0.npy' % name), spectrum_x_0, allow_pickle=True)
    np.save((config.HOME + 'outputs/%s_spectrum_x_1.npy' % name), spectrum_x_1, allow_pickle=True)
    np.save((config.HOME + 'outputs/%s_spectrum_x_2.npy' % name), spectrum_x_2, allow_pickle=True)
    np.save((config.HOME + 'outputs/%s_spectrum_y.npy' % name), spectrum_y, allow_pickle=True)
    np.save((config.HOME + 'outputs/%s_x_array.npy' % name), x_array, allow_pickle=True)

    return master_sz_0, master_sz_1, master_sz_2, spectrum_x_0, spectrum_x_1, spectrum_x_2, master_yt, ys, ts, realx_array, master_y, x_array, master_x_0, master_x_1, master_x_2, spectrum_y


def run_szpack(maps, sz_amp, band, isim, y_c, t_c, terp_x, terp_y, xaxis):
    nu = 300000.0 / clus_get_lambdas(maps[band]['band'])
    dI, thisx, JofXout, xout, errmsg = clus_get_relsz(isim, nu, band, y=y_c, te=t_c)
    sz_amp[band] = dI
    terp_x[band] = thisx
    terp_y[band] = JofXout
    xaxis[band] = xout


if __name__ == '__main__':
    sz_grid, input_yt, ys, ts = clus_szgf(0, 'rxj1347', 10, 0.1)
