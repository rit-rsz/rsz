################################################################################
# NAME : xid_test.py
# DATE STARTED : June 25, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE : Test for interfacing with XID
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import sys
sys.path.append('/home/vaughan/XID_plus/')
from xidplus import moc_routines
import xidplus
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
        data = hdul[1].data #convert to mJy, do we need to do this conversion? #this is map of data
        data2 = hdul[4].data  #convert to mJy #this is a noise map I guess? could be error or mask not sure.
        hdul.close()
        headers.append(header)
        datas.append(data)
        headers2.append(header2)
        datas2.append(data2)

    #getting a catalog thingy... nvm we ignore this.
    catdir = '/data/mercado/SPIRE/catalogs/'
    catfilename = 'srcext_rxj1347_1.00000.sav'
    inra = np.array([40.037706,       40.044160,       39.865923,       39.985885,       39.984657,       39.949583,       39.857725,       40.043405,
                     39.839410,       39.957867,       40.001735,       39.979558,       39.989626,       40.056271,       39.926356,       40.044352,
                     40.030643,       39.874660,       39.973949,       39.892760,       40.056794,       39.854375,       40.038465,       39.928479,
                     39.886157,       39.868784,       39.887424,       39.918426,       39.927220,       39.917336,       40.040241,       39.964287,
                     39.947827,       39.922271,       39.897630,       40.067808,       39.875642,       39.951659,       39.948114,       39.983048,
                     40.040930,       39.988071,       39.884689,       40.014071,       40.066620,       39.972957,       39.906640,       39.887340,
                     39.954942,       40.006865,       39.987952,       39.993509,       39.915210,       39.952706,       40.011619,       40.012867])
    indec = np.array([ -1.4987981,      -1.6229950,      -1.4773380,      -1.5524178,      -1.4239922,      -1.6170079,      -1.5154460,      -1.4291278,
                       -1.5406166,      -1.4487978,      -1.5868947,      -1.6795853,      -1.6706081,      -1.4690374,      -1.4853308,      -1.4536954,
                       -1.6434536,      -1.5542393,      -1.6120087,      -1.5311672,      -1.6558877,      -1.5971845,      -1.4270781,      -1.6189467,
                       -1.4691272,      -1.5430776,      -1.5772155,      -1.6189547,      -1.5423308,      -1.6353589,      -1.4396253,      -1.6401918,
                       -1.5588535,      -1.4962362,      -1.5739920,      -1.5381453,      -1.5141986,      -1.4838425,      -1.6755594,      -1.4956657,
                       -1.6309837,      -1.6567645,      -1.5454269,      -1.6270252,      -1.5774064,      -1.5505334,      -1.5525610,      -1.5383441,
                       -1.5808065,      -1.6343497,      -1.5622154,      -1.5681116,      -1.6522783,      -1.5686995,      -1.5949516,      -1.5732228])

    # c = SkyCoord(ra=[100]*u.degree, dec=[10]*u.degree) #i don't know if we need this or not
    # moc = pymoc.util.catalog.catalog_to_moc(c, 100, 15) #i don't know about the radius or resolution for the moc.
    # #this code seems to create a moc centered on the position given from the instance of the skycoords class.
    priors = []
    prf_size = np.array([18.15,25.15,36.3])
    prfs = []
    pixsize = np.array([maps[0]['pixsize'], maps[1]['pixsize'], maps[2]['pixsize']])
    pinds = []
    for i in range(len(files)):
        prior = xidplus.prior(datas[i], datas2[i], headers2[i], headers[i])
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
        print(priors[i].nsrc)
        priors[i].get_pointing_matrix()#calculate the pointing get_pointing_matrix
        #get_pointing_matrix(bkg=True) if bkg = True bkg is fitted to all pixels, else only fitted where prior souces contribute.
        print(priors[i].amat_data)
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
        print(key, fit.data[key])






if __name__ == '__main__':
    maps, err = get_data('a0370', verbose=1)
    xid_test(maps)
