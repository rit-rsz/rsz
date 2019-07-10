################################################################################
# NAME : get_xid.py
# DATE STARTED : June 21, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE : a script to interface with the XID program
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import numpy as np
import sys
sys.path.append('../utilities')
from get_spire_beam_fwhm import *
import matplotlib as plt
import config
from astropy.io import fits
from xidplus import moc_routines
import xidplus
import xidplus.catalogue as cat
import sys
sys.path.append('/home/vaughan/XID_plus/')
from xidplus import moc_routines
import xidplus
from scipy.io import readsav
sys.path.append('/home/vaughan/rsz/source_handling')
from get_data import *
import numpy as np
import pymoc
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.convolution import Gaussian2DKernel
from xidplus.stan_fit import SPIRE
from astropy.wcs import WCS as wcs
from astropy.wcs.utils import skycoord_to_pixel
import json
from astropy.coordinates import SkyCoord
from astropy import units as u


def get_xid(maps, cats, savemap=0, simmap=0, verbose=1, confusionerrors=1):
    err = False
    if simmap > 0:
        thresh = 3.0
        mf = 1.0
    else:
        thresh = 3.0
        mf = 1.0 #don't really know what the point of this if else was but ok.
    mJy2Jy = 1000.0 / mf
    catdir = config.CLUSDATA # this is a placeholder for now...
    catfile = config.CLUSDATA + 'placeholder' #this doesn't actually seem to do anything in the xid code,
    # but you need something for this for some reason.

    #this is dumb idk why I wrote this...
    # print('Retrieving data from cats')
    # inra = []
    # indec = []
    # for i in range(len(cats)):
    #     for j in range(len(cats[i]['ra'])):
    #         if cats[i]['ra'][j] not in inra and cats[i]['dec'][j] not in indec:
    #             inra.append(cats[i]['ra'][j])
    #             indec.append(cats[i]['dec'][j])
    # print('Done retrieving data from cats')

    inra = cats['ra']
    indec = cats['dec']


    #initializing data containers.
    pinds = []
    files = []
    primary_hdus = []
    noise_maps = []
    data_maps = []
    headers = []
    pixsizes = []
    prf_sizes = []
    priors = []
    prfs = []
    for i in range(len(maps)):
        #getting data from the fits files
        files.append(maps[i]['file'])
        hdul = fits.open(files[i])
        headers.append(hdul[1].header)
        primary_hdus.append(hdul[0].header)
        data_maps.append(hdul[1].data)
        noise_maps.append(hdul[4].data) #don't know if this is the error image or the mask image.
        pixsizes.append(maps[i]['pixsize'])
        prf_sizes.append(get_spire_beam_fwhm(maps[i]['band']))
        pinds.append(np.arange(0,101,1) * 1.0 / pixsizes[i]) #maybe this value needs to change?

        #setting up priors
        prior = xidplus.prior(data_maps[i], noise_maps[i], primary_hdus[i], headers[i])
        prior.prior_cat(inra, indec, catfile)
        prior.prior_bkg(-5.0, 5)

        #setting up prfs.
        prf = Gaussian2DKernel(prf_sizes[i] / 2.355, x_size=101, y_size = 101) #maybe x_size and y_size need to change.
        prf.normalize(mode='peak')
        prfs.append(prf.array)

        #appending prf to prior and setting point matrix
        prior.set_prf(prfs[i], pinds[i], pinds[i]) #prfs, xpinds, ypinds
        prior.get_pointing_matrix()
        prior.upper_lim_map()

        #appending prior to priors list.
        priors.append(prior)

    fit = SPIRE.all_bands(priors[0], priors[1], priors[2], iter=100) #number of iterations should be at least 100 just set lower for testing.
    posterior = xidplus.posterior_stan(fit,[priors[0],priors[1],priors[2]])

    figs, fig = xidplus.plots.plot_Bayes_pval_map(priors, posterior)
    # print(type(figs)) #figs is list.
    # print(figs) #fig is matplotlib.figure.figure object.
    # print(type(fig))
    cols = ['PSW', 'PMW', 'PLW']
    counter = 0
    for figure in figs:
        figure.save('xid_%s.png' %(cols[counter]))
        counter += 1

    # plt.imshow(figs)

    spire_cat = cat.create_SPIRE_cat(posterior, priors[0], priors[1], priors[2])

    xid_data = spire_cat[1].data
    xid = []

    #in units of mJy for fluxes and degrees for RA/DEC
    xid1 = {#'sid' : priors[0].ID,
            'band' : 'PSW',
            'sra' : xid_data.field('RA'),
            'sdec' : xid_data.field('DEC'),
            'sflux' : xid_data.field('F_SPIRE_250'),
            'serr' : xid_data.field('FErr_SPIRE_250_u'), #there was also FErr_SPIRE_250_l don't know which to use.
            'pflux' : xid_data.field('F_SPIRE_250'),
            'perr' : xid_data.field('FErr_SPIRE_250_u'),
            'x' : None,
            'y' : None,
            'model' : None,
            'mf' : mf} #idk if perr and pflux is right there may be a conversion needed for pflux.
            #in mikes code it has pflux = output from xid / mJy2Jy.
    xid2 = {#'sid' : priors[1].ID,
            'band' : 'PMW',
            'sra' : xid_data.field('RA'),
            'sdec' : xid_data.field('DEC'),
            'sflux' : xid_data.field('F_SPIRE_350'),
            'serr' : xid_data.field('FErr_SPIRE_350_u'),
            'pflux' : xid_data.field('F_SPIRE_350'),
            'perr' : xid_data.field('FErr_SPIRE_350_u'),
            'x' : None,
            'y' : None,
            'model' : None,
            'mf' : mf}

    xid3 = {#'sid' : priors[2].ID,
            'band' : 'PLW',
            'sra' : xid_data.field('RA'),
            'sdec' : xid_data.field('DEC'),
            'sflux' : xid_data.field('F_SPIRE_500'),
            'serr' : xid_data.field('FErr_SPIRE_500_u'),
            'pflux' : xid_data.field('F_SPIRE_500'),
            'perr' : xid_data.field('FErr_SPIRE_500_u'),
            'x' : None,
            'y' : None,
            'model' : None,
            'mf' : mf}

    #there was another term in the dictionary sflux, pflux and sflux looked like maybe the same thing, but I'm not sure.
    #I left it out so if there are issues with that then it is because that is gone.
    xid.append(xid1)
    xid.append(xid2)
    xid.append(xid3)

    # only look at data with a flux lower than 0.0
    for i in range(len(xid)):
        whpl = []
        for j in range(xid[i]['pflux'].shape[0]):
            if xid[i]['pflux'][j] >= 0.0:
                whpl.append(j)
        whpl = np.array(whpl)

        xid[i]['sra'] = xid[i]['sra'][whpl]
        xid[i]['sdec'] = xid[i]['sdec'][whpl]
        xid[i]['x'] = xid[i]['x'][whpl]
        xid[i]['y'] = xid[i]['y'][whpl]
        xid[i]['sflux'] = xid[i]['sflux'][whpl]
        xid[i]['serr'] = xid[i]['serr'][whpl]

          #initializing w class.

    # w = wcs(spire_cat[1].header)

    #converting data to a list instead of numpy array so that we can save it in a .json file.
    for i in range(len(xid)):
        ra = xid[i]['sra'] * u.deg
        dec = xid[i]['sdec'] * u.deg
        c = SkyCoord(ra, dec)
        #initializing w class.
        hdul = fits.open(maps[i]['file'])
        w = wcs(hdul[1].header)
        #converting ra/dec to pixel coords.
        px, py = skycoord_to_pixel(c, w)
        xid[i]['x'] = px
        xid[i]['y'] = py
        xid[i]['sra'] = xid[i]['sra'].tolist()
        xid[i]['sdec'] = xid[i]['sdec'].tolist()
        xid[i]['sflux'] = xid[i]['sflux'].tolist()
        xid[i]['serr'] = xid[i]['serr'].tolist()
        xid[i]['pflux'] = xid[i]['pflux'].tolist()
        xid[i]['perr'] = xid[i]['perr'].tolist()
        xid[i]['x'] = xid[i]['x'].tolist()
        xid[i]['y'] = xid[i]['y'].tolist()

        #saving to json file for further analysis.
        with open('xid_%s.json' %(xid[i]['band']), 'w') as f: #code for saving output to a file.
            json.dump(xid[i], f)

    #model = image_model(x,y, sflux, maps[i]['astr']['NAXIS'][0], maps[i]['astr']['NAXIS'][1],
    #maps[i]['psf'])
    #need to finish converting model over to python.











        if savemap:
            outfile = config.CLUSSBOX + 'clus_get_xid_model_' + maps[i]['band'] + '.fit'
            writefits(outfile, data=model, header_dict=maps[i]['shead'])



    return xid, err
