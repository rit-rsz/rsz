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
        header = hdul[1].header #map hdu
        header2 = hdul[0].header #primary hdu
        data = hdul[1].data * 1*10**3 #convert to mJy, do we need to do this conversion? #this is map of data
        data2 = hdul[2].data * 1*10**3 #convert to mJy #this is a noise map I guess? could be error or mask not sure.
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

    c = SkyCoord(ra=[100]*u.degree, dec=[10]*u.degree) #i don't know if we need this or not
    moc = pymoc.util.catalog.catalog_to_moc(c, 100, 15) #i don't know about the radius or resolution for the moc.
    #this code seems to create a moc centered on the position given from the instance of the skycoords class.
    priors = []
    prf_size = np.array([18.15,25.15,36.3])
    prfs = []
    pixsize = np.array([maps[0]['pixsize'], maps[1]['pixsize'], maps[2]['pixsize']])
    pinds = []
    for i in range(len(files)):
        prior = xidplus.prior(datas[i], datas2[i], headers[i], headers2[i])
        prior.prior_cat(inra, indec, catfilename)
        #not clear what this does?
        #prior_cat(ra, dec, prior_cat_file, flux_lower, flux_upper, id, moc, z_median, z_sig)
        #ra - right ascension of sources, dec - declination of sources, prior_cat_file - catfile
        #flux_lower - lower limit of flux for each source flux_upper - upper limit of flux for each source
        #ID - help id for each source, moc - multi-order-coverage map z_median - mean of redshift pdf
        #z_sig - sigma of redshift pdf
        prior.prior_bkg(-5.0, 5) #sets mean and sigma. should these be different values?
        #prior_bkg(mu, sigma) mu = mean, sigma = standard deviation
        priors.append(prior)
        prf = Gaussian2DKernel(prf_size[i] / 2.355, x_size=101, y_size=101)
        #create a gaussian based off of the FWHM of each of the different beams
        prf.normalize(mode='peak')
        #normalize the gaussian kernel for each beam? maybe this isn't necessary, not sure which normalize method? can't find documentation for it...
        prfs.append(prf.array)
        pind = np.arange(0,101,1) * 1.0 / pixsize[i] #get scale in terms of pixel map
        pinds.append(pind)

    for i in range(len(priors)):
        priors[i].set_prf(prfs[i],pinds[i],pinds[i]) #add prf array and corresponding x and y scales.
        #set_prf(prf, pindx, pindy) prf - nxn array and center or prf is center of array
        #pindx - n array pixel scale of prf array #pindy - n array pixel scale of prf array.
        #whats the difference between pindx and pindy?
        priors[i].get_pointing_matrix()#calculate the pointing get_pointing_matrix
        #get_pointing_matrix(bkg=True) if bkg = True bkg is fitted to all pixels, else only fitted where prior souces contribute.
        priors[i].upper_lim_map() #updates the flux upper limit to abs(bkg) + 2*sigma_bkg + max(D) where max(D) is the maximum value of pixels the source contributes to.

    fit = SPIRE.all_bands(priors[0], priors[1], priors[2], iter=1000)
    #fits the three spire bands.
    #SPIRE.all_bands(spire1, spire2, spire3, chains, iter)
    #spire1 - spire image 1
    #spire2 - spire image 2
    #spire3 - spire image 3
    #chains - the number of chains
    #iter - # the number of iterations
    for key in fit.data.keys():
        print(key)






if __name__ == '__main__':
    maps, err = get_data('rxj1347', verbose=1)
    xid_test(maps)
