import numpy as np
import config

def get_tfs(clusname):
    tf_file = config.clus_path + clusname + '_itermap_tf.fits' #missing !CLUSSBOX don't know
    tf_maps, err = clus_get_data(clusname,manpath=tf_file)
    if err:
        print('error when running clus_get_data: ', err)
    else:
        ncols = len(tf_maps)
        for i in range(ncols):
            tf_maps[i].xclean = tf_maps[i].signal
            tf_maps[i].error = np.tile(1, tf_maps[i].astr.naxis) #don't know what .astr.naxis is yet
            tf_maps[i].calfac = 1./tf_maps.JY2MJy
    return tf_maps
