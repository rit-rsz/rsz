
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
sys.path.append('../multiband_pcat')
from pcat_spire import lion

def clus_pcat_setup(maps,params,err=None):

    print('working')
    # ob = lion(raw_counts=True, auto_resize=True, visual=True)
    ob = lion(map_object=maps,band1=1, auto_resize=True, visual=True, return_median_model=True, make_post_plots=True, nsamp=1, residual_samples=1)
    x = ob.main()
    # Have a commented list of all the settings here for testing usage.
    #Need to figure out the order of all the args

    exit()
    return maps,err

def clus_pcat_setup_python2():

    ob = lion(band1=1, auto_resize=True, visual=True, make_post_plots=True, nsamp=1000, residual_samples=200)
    ob.main()
    exit()

if __name__ == '__main__':
    clus_pcat_setup_python2()
