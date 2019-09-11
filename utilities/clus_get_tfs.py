################################################################################
# NAME : clus_get_tfs
# DATE STARTED : June 20, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE : create a python counterpart to clus_get_tfs.pro
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################


import numpy as np
import config
import sys
sys.path.append('source_handling')
from clus_get_data import *

def clus_get_tfs(clusname):
    err = False
    tf_dir = config.CLUSSBOX
    tf_id =  clusname + '_itermap_tf'
    tf_maps, err = clus_get_data(clusname,manpath=tf_dir, manidentifier=tf_id)
    if err:
        print('error when running clus_get_data: ', err)
    else:
        ncols = len(tf_maps)
        for i in range(ncols):
            tf_maps[i]['xclean'] = tf_maps[i]['signal']
            tf_maps[i]['error'] = np.tile(1, tf_maps[i]['astr']['NAXIS']) #don't know what .astr.naxis is yet
            tf_maps[i]['calfac'] = 1./tf_maps['JY2MJy']
    return tf_maps, err
