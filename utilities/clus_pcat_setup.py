
################################################################################
# NAME : add_sziso.py
# DATE STARTED : September 18, 2019
# AUTHORS : Dale Mercado
# PURPOSE : Used to set up the inputs for Pcat to called and fed our signal maps
# EXPLANATION : There are a few things that need to be passed into pcat inorder for
#               it to work as intended. To do this we need to hand it the signal map,
#               the error map as well as some other setting parameters.
#               the issue is that we currently would have to pass each in, unless
#               we can get it to accept the map object and parse it in pcat
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
#           pcat
#           err
#           These will probably change and are based on how xid was original set up
# REVISION HISTORY :
################################################################################

import sys
import config
sys.path.append(config.HOME + 'multiband_pcat')
from pcat_spire import lion
import numpy as np

def clus_pcat_setup(maps,params,err=None):

    # ob = lion(raw_counts=True, auto_resize=True, visual=True)
    ob = lion(band0=0, map_object=maps, auto_resize=True, make_post_plots=False, nsamp=10, residual_samples=10)
    x = ob.main()

    return x,err

if __name__ == '__main__':
    clus_pcat_setup_python2()
