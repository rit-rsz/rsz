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

import matplotlib.pyplot as plt
from xidplus import moc_routines
import xidplus
import xidplus.catalogue as cat
import sys
from scipy.io import readsav
sys.path.append('/home/vaughan/rsz/source_handling')
from get_data import *
from astropy.io import fits
import numpy as np
import pymoc
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.convolution import Gaussian2DKernel
from xidplus.stan_fit import SPIRE
from astropy.wcs import WCS as wcs

def xid_test(maps):
    mJy2Jy = 1000.0
    files = []
    # for i in range(len(maps)):
    #     files.append(maps[i]['file'])
    file = '/data/mercado/SPIRE/bethermin_sims/a0370/a0370_PSW_sim0200.fits'
    files.append(file)
    file = '/data/mercado/SPIRE/bethermin_sims/a0370/a0370_PMW_sim0200.fits'
    files.append(file)
    file = '/data/mercado/SPIRE/bethermin_sims/a0370/a0370_PLW_sim0200.fits'
    files.append(file)
    headers = []
    datas = []
    datas2 = []
    headers2 = []
    for file in files:
        hdul = fits.open(file)
        header = hdul[1].header #map hdu
        header2 = hdul[0].header #primary hdu
        data = hdul[1].data #convert to mJy, do we need to do this conversion? #this is map of data
        data2 = hdul[3].data  #convert to mJy #this is a noise map I guess? could be error or mask not sure.
        hdul.close()
        headers.append(header)
        datas.append(data)
        headers2.append(header2)
        datas2.append(data2)
        hdul.close()

    pixsizes = []
    # for i in range(len(maps)):
    #     pixsizes.append(maps[i]['pixsize'])


    for header in headers:
        pixsize = 3600 * \
                  np.mean([abs(header['CD1_1']+header['CD2_1']), \
                         abs(header['CD2_1'] + header['CD2_2'])])
        pixsizes.append(pixsize)
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
    pixsize = np.array([pixsizes[0], pixsizes[1], pixsizes[2]])
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
        priors[i].set_prf(prfs[i],pinds[i],pinds[i]) #add prf a    pixel_coords = []
        #set_prf(prf, pindx, pindy) prf - nxn array and center or prf is center of array
        #pindx - n array pixel scale of prf array #pindy - n array pixel scale of prf array.
        priors[i].get_pointing_matrix()#calculate the pointing get_pointing_matrix
        #get_pointing_matrix(bkg=True) if bkg = True bkg is fitted to all pixels, else only fitted where prior souces contribute.
        #priors[i].upper_lim_map() #updates the flux upper limit to abs(bkg) + 2*sigma_bkg + max(D) where max(D) is the maximum value of pixels the source contributes to.

    fit = SPIRE.all_bands(priors[0], priors[1], priors[2], iter=100)
    posterior = xidplus.posterior_stan(fit,[priors[0],priors[1],priors[2]])

    spire_cat = cat.create_SPIRE_cat(posterior, priors[0], priors[1], priors[2])
    print(fit)
    print(spire_cat.info())

    # w = wcs('test1.fits')

    # for key in spire_cat[1].header.keys():
    #     print(key)
    ra = spire_cat[1].data.field('RA')
    dec = spire_cat[1].data.field('DEC')

    x, y = w.all_pix2world(ra, dec, 1.0)
    print(len(inra))
    print(len(indec))
    print(len(ra))
    print(len(dec))

    plt.scatter(x,y)
    plt.show()

    print(posterior.sra)
    # plt.plot(figs)
    # print(fit.summary('src_f')['summary'])
    #fits the three spire bands.
    #SPIRE.all_bands(spire1, spire2, spire3, chains, iter)
    #spire1 - spire image 1
    #spire2 - spire image 2
    #spire3 - spire image 3
    #chains - the number of chains
    #iter - # the number of i    hdul = fits.open('/data/mercado/SPIRE/hermes_clusters/a0370_PLW_fr_1.fits')
    # print(hdul.info())
    # print(fit)
    # print(fit.summary())
    pixels = fit.data['Row_psw'] #the pixel contributing to the source number, I don't get how this is related to the actual pixel number.
    sources = fit.data['Col_plw'] #what source number we are on
    val = fit.data['Val_plw']

    # for i in range(len(sources)):
    #     indexes = np.where(sources == i)
    #     star_pixels = pixels[indexes]
    #
    #
    #
    # plt.scatter(px,py)
    # plt.show()
    #





if __name__ == '__main__':
    maps, err = get_data('rxj1347')
    xid_test(maps)
