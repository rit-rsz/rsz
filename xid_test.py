from xidplus import moc_routines
import xidplus
import sys
from scipy.io import readsav
sys.path.append('source_handling')
from get_data import *
from astropy.io import fits
import numpy as np
import pymoc
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.convolution import Gaussian2DKernel
from xidplus.stan_fit import SPIRE

def xid_test(maps):
    length = len(maps)
    files = []
    for i in range(length):
        files.append(maps[i]['file'])
    headers = []
    datas = []
    datas2 = []
    headers2 = []
    for file in files:
        hdul = fits.open(file)
        header = hdul[1].header
        header2 = hdul[0].header
        data = hdul[1].data * 1*10**3 #convert to mJy
        data2 = hdul[4].data * 1*10**3 #convert to mJy
        hdul.close()
        headers.append(header)
        datas.append(data)
        headers2.append(header2)
        datas2.append(data2)

    #getting a catalog thingy... nvm we ignore this.
    catdir = '/data/mercado/SPIRE/catalogs/'
    catfilename = 'srcext_rxj1347_1.00000.sav'
    inra = np.arange(90)
    indec = np.arange(90)

    # c = SkyCoord(ra=[100]*u.degree, dec=[10]*u.degree) #can maybe have skycoords a list of RAs/ DECs
    # moc = pymoc.util.catalog.catalog_to_moc(c, 100, 15) #don't think we need a moc but not sure...

    priors = []
    prf_size = np.array([18.15,25.15,36.3])
    prfs = []
    pixsize = np.array([maps[0]['pixsize'], maps[1]['pixsize'], maps[2]['pixsize']])
    pinds = []
    for i in range(len(files)):
        prior = xidplus.prior(datas[i], datas2[i], headers[i], headers2[i])
        prior.prior_cat(inra, indec, catfilename)
        prior.prior_bkg(-5.0, 5)
        priors.append(prior)
        prf = Gaussian2DKernel(prf_size[i] / 2.355, x_size=101, y_size=101)
        prf.normalize(mode='peak')
        prfs.append(prf)
        pind = np.arange(0,101,1) * 1.0 / pixsize[i]
        pinds.append(pind)

    for i in range(len(priors)):
        priors[i].set_prf(prfs[i],pinds[i],pinds[i])
        priors[i].get_pointing_matrix()
        priors[i].upper_lim_map()

    fit = SPIRE.all_bands(priors[0], priors[1], priors[2], iter=1000)
    for key in fit.data.keys():
        print(fit.data[key])






if __name__ == '__main__':
    maps, err = get_data('rxj1347', verbose=1)
    xid_test(maps)
