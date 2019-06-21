"""
Name: clus_get_tfs
Date Started: June 20, 2019
Authors: Benjamin Vaughan
Purpose: create a python counterpart to clus_get_tfs.pro
Explanation:
Calling Sequence:
Inputs: clusname - name of the cluster
Outputs: tf_maps -
Revision History:
"""


import numpy as np
import config
import sys
sys.path.append('source_handling')
from get_data import *

def get_tfs(clusname):
    err = False
    tf_file = config.CLUSBOS
    tf_identifier = '_itermap_tf' #missing !CLUSSBOX don't know
    tf_maps, err = get_data(clusname,manpath=tf_file, manidentifier=tf_identifier)
    print(type(tf_maps[0]))
    if err:
        print('error when running clus_get_data: ' + err)
        return None, err
    else:
        ncols = len(tf_maps)
        for i in range(ncols):
            tf_maps[i].xclean = tf_maps[i].signal
            tf_maps[i].error = np.tile(1, tf_maps[i].astr.naxis) #don't know what .astr.naxis is yet
            tf_maps[i].calfac = 1./tf_maps.JY2MJy
    return tf_maps, err
