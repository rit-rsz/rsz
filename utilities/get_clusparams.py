
################################################################################
# NAME : get_clusparams.py
# DATE STARTED : June 11, 2019
# AUTHORS : Victoria Butler & Dale Mercado
# PURPOSE : Fetches the cluster lookup table
#           This should just be some of the header material
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import scipy.io
import numpy as np
# from config import * #(this line will give me access to all directory variables)
import matplotlib.pyplot as plt
import math
from astropy.io.ascii import read
from astropy.table import Column
# tab = read('spectrum.csv')
import csv
from collections import defaultdict


def get_clus_params(clusname,VERBOSE=verbose,SUCCESS=cgp_success,\
                           ERRMSG=cgp_errmsg):
    # verbose = 1 if not verbose else verbose

    clusname_in = str('rxj1347')
    tab = read('/data/mercado/SPIRE/' + 'lookup/cluster_master_list_150327.csv')
    # ,header_start = 0)
    # print(tab[';Cluster Name'][1])
    x = tab[';Cluster Name']

    # Header information
    fielddesc = ['clusname: clus descriptor',\
                   'longname: official NED field name',\
                   'program: program cluster belongs to (HLS or HerMES)',\
                   'z: redshift of source',\
                   'nomrah: nominal RA in hh:mm:ss',\
                   'nomrad: nominal RA in deg',\
                   'fidrah: fiducial SZ RA in hh:mm:ss',\
                   'fidrad: fiducial SZ RA in deg',\
                   'nomdeh: nominal Dec in dd:am,as',\
                   'nomded: nominal Dec in deg',\
                   'fiddeh: fiducial SZ Dec in dd:am:as',\
                   'fidded: fiducial SZ Dec in deg',\
                   'Te: temperature of ICM',\
                   'rc: isothermal beta model core radius',\
                   'beta: isothermal beta model slope',\
                   'yzero: central y parameter',\
                   'yup: upper uncertainty on yzero',\
                   'ydn: lower uncertainty on yzero',\
                   'cirrus: 1 = cirrus contamination present']

    # print(type(clusname_in))
    l = 57
    n = 0
    for i in range(0,l):
        y = str(x[i])
        # print(y)
        # print(clusname_in)
        # matchclus = np.equal(y, clusname_in)
        if y == clusname_in :
            n = i
            # print(n)
#           This is my work around to finding the index where clusname_in and clusname from the
#             input match up
    if n !=1:
        errmsg = 'Number of cluster name matches is not 1'
        # if verbose:
        #     print(errmsg)
        #     return(0)


    clusparams = np.array([fielddesc, tab[';Cluster Name'][n] ,\
    tab['Full Name'][n] ,\
    tab['Program'][n] ,\
    tab['z'][n] ,\
    tab['Nom. RA'][n] ,\
    tab['ND RA'][n] ,\
    tab['Fid. RA'][n] ,\
    tab['FD RA'][n] ,\
    tab['Nom Dec.'][n] ,\
    tab['ND Dec.'][n] ,\
    tab['Fid Dec.'][n] ,\
    tab['FD Dec.'][n] ,\
    tab['Te'][n] ,\
    tab['r_c'][n] ,\
    tab['beta'][n] ,\
    tab['y_0'][n] ,\
    tab['yup'][n] ,\
    tab['ydn'][n] ,\
    tab['cirrus'][n]])
    # dtype below is trying to recreate clusparams.clusname style
    # dtype=[('clusname'),('longname'),('program'),('z'),('nomrah'),\
    #         ('nomrad'),('fidrah'),('fidrad'),('nomdeh'),('nomded'),\
    #         ('fiddeh'),('fidded'),('Te'),('rc,beta'),('yzero'),('yup'),('ydn'),('cirrus')])
#   This currently is a gross solution to appending everything to get_clus_params
#   There must be cleaner way than this but this will do for the quick and dirty mean time.
#   One thing that this is missing is any header information.
#   In later scrpts the pipeline calls for instance params.fidrad and that is not a
#   feature that this current setup has which is important

    print(clusparams[2])

    #now need to set clusparams to be handed back to the script
    #Need to figure out what kind of variable to do this as

    return(clusparams)


if __name__ == '__main__':
    get_clus_params()
