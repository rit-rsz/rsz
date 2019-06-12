
################################################################################
# NAME : catsrc.py
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


def get_clus_params():

    # verbose = 1 if not verbose else verbose


    tab = read('/data/mercado/SPIRE/' + 'lookup/cluster_master_list_150327.csv')
    # ,header_start = 0)
    print(tab[';Cluster Name'])
    col = Column(name=tab[';Cluster Name'], dtype=str, length = 56)
    print(col)


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
    print(type(tab[';Cluster Name']))

    # matchclus = np.equal(tab[';Cluster Name'],clusname_in)
    # print(matchclus)
    # whpl = where(matchclus,count)

    # if count != 1:
    #     errmsg = 'Number of cluster name matches is not 1'
    #     if verbose:
    #         print errmsg

    #now need to set clusparams to be handed back to the script


if __name__ == '__main__':
    get_clus_params()
