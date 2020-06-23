################################################################################
# NAME : clus_sz_grids
# DATE STARTED : June 17, 2020
# AUTHORS :  Benjamin Vaughan
# PURPOSE : This routine is meant to create look up grids for the SZe which we
# use to do our spectral fitting.
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import argparse
import numpy as np
import sys
sys.path.append('/home/bjv7945/rsz/sz/')
sys.path.append('/home/bjv7945/rsz/utillities/')
from clus_get_relsz import *
from clus_get_lambdas import *


def thermal(v, I, y):

    func1 = (((v**4)*np.exp(v))/(np.exp(v)-1)**2) * ((v*((np.exp(v) + 1)/(np.exp(v) - 1)))-4)
    return func1*(I*y)


def create_grid(name, grid_hw):
    '''
    Inputs - clusname (str): the name of the cluster
           - grid_hw (int): the height and width of the grid (it is assumed to be a square.)
    Output: - a grid containing sz values
            - a grid containing spectra at those values
            - an array of ys values
            - an array of ts values
    '''
    bands = ['BOLOCAM', 'PSW', 'PMW', 'PLW']

    T_0     = 2.725                   #CMB temperature, K
    k_B     = 1.3806503e-23           #Boltzmann constant, J/K
    h       = 6.626068e-34            #Planck constant, J s
    c       = 2.99792458e8            #m/s
    GHztoHz = 1.0e9                   #GHz -> Hz

    #later on will want to not have the intervals hard coded into the script.
    y_in = np.linspace(0, 12e-4, grid_hw)
    t_in = np.linspace(0, 25, grid_hw)

    #this is some math to calculate the cannonical effect when Te is set to zero.
    I = 2*((k_B*T_0)**3/(h*c)**2) * 10**(26) * 10**(-6)
    v = (h/(k_B * T_0))*np.linspace(1*10**(9),1700*10**9,100)

    print('Grid size is (%s, %s) with %s points' % (len(y_in), len(t_in), len(y_in)**2))
    #create empty data structure
    bolo_sz_grid = np.zeros((len(y_in), len(t_in)))
    psw_sz_grid = np.zeros((len(y_in), len(t_in)))
    pmw_sz_grid = np.zeros((len(y_in), len(t_in)))
    plw_sz_grid = np.zeros((len(y_in), len(t_in)))
    spec_y  = np.zeros((len(y_in) , len(t_in), 100))
    spec_x  = np.zeros((len(y_in) , len(t_in), 100))
    y_t_grid = np.zeros((len(y_in), len(t_in), 2))
    #iterate through the y and T values to create a grid of the SZe
    #in this loop we only use the bolocam runs to save the spectrum to reduce storage cost.
    for i in range(len(y_in)):
        for j in range(len(t_in)):
            if t_in[j] == 0: #check for Te = 0, in this case we want to use the cannonical effect.
                spec = thermal(v, I, y_in[i])
                spec_y[i, j, :] = spec
                spec_x[i, j, :] = v
                wavelen = [clus_get_lambdas(bandi) for bandi in bands]
                v_band = [h * waveleni / (k_B * T_0) for waveleni in wavelen]
                dI = [thermal(v_bandi, I, y_in[i]) for v_bandi in v_band]
                bolo_sz_grid[i, j] = dI[0]
                psw_sz_grid[i, j] = dI[1]
                pmw_sz_grid[i, j] = dI[2]
                plw_sz_grid[i, j] = dI[3]
                y_t_grid[i , j, 0] = y_in[i]
                y_t_grid[i , j, 1] = t_in[j]
            else:
                dI, x_pos, spectrum, xaxis = run_szpack('ALL', y_in[i], t_in[j])
                bolo_sz_grid[i, j] = dI[0]
                psw_sz_grid[i, j] = dI[1]
                pmw_sz_grid[i, j] = dI[2]
                plw_sz_grid[i, j] = dI[3]
                y_t_grid[i , j, 0] = y_in[i]
                y_t_grid[i , j, 1] = t_in[j]
                spec_y[i, j, :] = spectrum
                spec_x[i, j, :] = xaxis

    #save the output arrays.
    np.save(config.OUTPUT + 'sz_grids/%s_sz_grid_BOLO.npy' % (name), bolo_sz_grid, allow_pickle=True)
    np.save(config.OUTPUT + 'sz_grids/%s_sz_grid_PSW.npy' % (name), psw_sz_grid, allow_pickle=True)
    np.save(config.OUTPUT + 'sz_grids/%s_sz_grid_PMW.npy' % (name), pmw_sz_grid, allow_pickle=True)
    np.save(config.OUTPUT + 'sz_grids/%s_sz_grid_PLW.npy' % (name,), plw_sz_grid, allow_pickle=True)
    np.save(config.OUTPUT + 'sz_grids/%s_ys_grid.npy' % name, y_in, allow_pickle = True)
    np.save(config.OUTPUT + 'sz_grids/%s_ts_grid.npy' % name, t_in, allow_pickle = True)
    np.save(config.OUTPUT + 'sz_grids/%s_y_t_grid.npy' % name, y_t_grid, allow_pickle=True)

    np.save((config.OUTPUT + 'sz_grids/%s_spectrum_y.npy' % name), spec_y, allow_pickle=True)
    np.save((config.OUTPUT + 'sz_grids/%s_spectrum_x.npy' % name), spec_x, allow_pickle=True)

def run_szpack(band, y_c, t_c):
    '''
    This is a wrapper to run the SZpack script.
    Inputs - band (str) : the band we want to find the dI for, must be ('PSW', 'PMW', 'PLW', or 'BOLOCAM') 'ALL' returns all of these.
           - y_c (float): the input value for the compton y parameter
           - t_c (float): the input value for the Temperature
    Outputs - dI (float) : the intensity in Mjy / Sr of the SZe at the given band
            - thisx (float) : the xaxis coordinate of dI on the output spectra.
            - xout (float array) : an array of xpositions for the spectra (in units of dimensionless freq see Carlstrom et al. 2012)
            - JofXout (float array) : this is an array of the SZe spectra.
    '''
    #run SZpack.
    dI, thisx, JofXout, xout, errmsg = clus_get_relsz(0, 0, band, y=y_c, te=t_c)
    return dI, thisx, JofXout, xout


if __name__ == '__main__':
    '''
    The following code is for use with the SPORC super computer so that we can parallize this process.
    '''
    print('Creating grids for likelihood analysis')
    parser = argparse.ArgumentParser()
    parser.add_argument("-run", help="This runs real or sim map through the pipeline for a single map on a specific cluster.", nargs=2, metavar=('clusname', 'grid_size'))
    args = parser.parse_args()
    if args.run:
        clusname  = args.run[0]
        gridding = args.run[1]
        create_grid(clusname, gridding)
    #create a grid of y and t values based off of the nominal value from the cluster
