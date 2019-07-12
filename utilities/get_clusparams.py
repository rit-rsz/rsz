################################################################################
# NAME : get_clusparams.py
# DATE STARTED : June 11, 2019
# AUTHORS : Dale Mercado
# PURPOSE : Fetches the cluster lookup param_data
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
from astropy.io.ascii import read
import config


def get_clus_params(clusname_in,
                    verbose = 1):#, su    return errmsg
    errmsg = False
#   Read in the csv file to work with
    param_data = read(config.CLUSDATA + 'lookup/cluster_master_list_150327.csv')
#   Since cluster name is an easy search we choose to set it early
#   Needs to be set like this because calling param_data[1] doesn't seem to work as intended
    param_name = param_data[';Cluster Name']

#   Header information
    fielddesc = ['clusname: clus descriptor',\
                   'longname: official NED field name',\
                   'program: program cluster belongs to (HLS or HerMES)',\
                   'z: redshift of source',\
                   'nomrah: nominal RA in hh:mm:ss',\
                   'nomrad: nominal RA in deg',\
                   'fidrah: fiducial SZ RA in hh:mm:ss', # the actual measured position of the SZ effect
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

#   This is the best python equivalent that I could conseve to counter idl's where function
    param_len = len(param_data)
    count = 0
    for i in range(param_len):
        place = str(param_name[i])

        if place == clusname_in :
            whpl = i
            count += 1
#           Where place, This will indicate the position of the cluster in the sav file,
#           this then gets passed onto the clusparams dictionary for further use


    if count !=1:
        errmsg = 'Number of cluster name matches is not 1'
        return None, errmsg

#   Creates a dictionary for the selected galaxy and returns it back to catsource
    clusparams = {'fielddesc':fielddesc,\
    'clusname':param_data[';Cluster Name'][whpl] ,\
    'longname':param_data['Full Name'][whpl] ,\
    'program':param_data['Program'][whpl] ,\
    'z':param_data['z'][whpl] ,\
    'nomrah':param_data['Nom. RA'][whpl] ,\
    'nomrad':param_data['ND RA'][whpl] ,\
    'fidrah':param_data['Fid. RA'][whpl] ,\
    'fidrad':param_data['FD RA'][whpl] ,\
    'nomdeh':param_data['Nom Dec.'][whpl] ,\
    'nomded':param_data['ND Dec.'][whpl] ,\
    'fiddeh':param_data['Fid Dec.'][whpl] ,\
    'fidded':param_data['FD Dec.'][whpl] ,\
    'Te':param_data['Te'][whpl] ,\
    'rc':param_data['r_c'][whpl] ,\
    'beta':param_data['beta'][whpl] ,\
    'yzero':param_data['y_0'][whpl] ,\
    'yup':param_data['yup'][whpl] ,\
    'ydin':param_data['ydn'][whpl] ,\
    'cirrus':param_data['cirrus'][whpl]}


#   print(clusparams['clusname'])

#   Hand off clusparams back to catsrc
    return clusparams, errmsg

#For debugging, used to call this file directly
if __name__ == '__main__':
    get_clus_params('rxj1347')
