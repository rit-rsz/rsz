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
sys.path.append('utilities')
from get_spire_beam_fwhm import *

def get_xid(maps, cat, savemap=0, simmap=0, verbose=1, confusionerrors=1):
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

    #data from cat
    inra = cat['ra']
    indec = cat['dec']

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
    for i in range(len(map)):
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
        priors.append(prior)\

    fit = SPIRE.all_bands(priors[0], priors[1], priors[2], iter=10) #number of iterations can be changed easily.

    xid = []

    xid1 = {'sid' : 1,
            'band' : 'PSW',
            'sra' : inra,
            'sdec' : indec,
            'sflux' : fit['Val_psw'],
            'serr' : fit['sigma_psw'],
            'pflux' : fit['Val_psw'],
            'perr' : fit['sigma_psw'],
            'x' : fit['Row_psw'],
            'y' : fit['Col_psw'],
            'model' : None,
            'mf' : mf} #idk if perr and pflux is right there may be a conversion needed for pflux.
            #in mikes code it has pflux = output from xid / mJy2Jy.
    xid2 = {'sid' : 2,
            'band' : 'PMW',
            'sra' : inra,
            'sdec' : indec,
            'sflux' : fit['Val_pmw'],
            'serr' : fit['sigma_pmw'],
            'pflux' : fit['Val_pmw'],
            'perr' : fit['sigma_pmw'],
            'x' : fit['Row_pmw'],
            'y' : fit['Col_pmw'],
            'model' : None,
            'mf' : mf}

    xid3 = {'sid' : 3,
            'band' : 'PLW',
            'sra' : inra,
            'sdec' : indec,
            'sflux' : fit['Val_plw'],
            'serr' : fit['sigma_plw'],
            'pflux' : fit['Val_plw'],
            'perr' : fit['sigma_plw'],
            'x' : fit['Row_plw'],
            'y' : fit['Col_plw'],
            'model' : None,
            'mf' : mf}

    #there was another term in the dictionary sflux, pflux and sflux looked like maybe the same thing, but I'm not sure.
    #I left it out so if there are issues with that then it is because that is gone.
    xid.append(xid1)
    xid.append(xid2)
    xid.append(xid3)
    #don't know if inra gets changed or the input is just spit back out, but in the python XID version there is
    #Row_psw and Col_psw for all bands, I think this refers to the pixel coordinates of sources within that band.

    #model = image_model(x,y, sflux, maps[i]['astr']['NAXIS'][0], maps[i]['astr']['NAXIS'][1],
    #maps[i]['psf'])
    #need to finish converting model over to python.


    for i in range(len(xid)):
        whpl = []
        for j in range(xid[i]['pflux'].shape[0]):
            for k in range(xid[i]['pflux'].shape[1]):
                if xid[i]['pflux'][j,k] >= 0.0:
                    whpl.append([j,k])
        whpl = np.array(whpl)

        xid[i]['sra'] = xid[i]['sra'][whpl]
        xid[i]['sdec'] = xid[i]['sdec'][whpl]
        xid[i]['x'] = xid[i]['x'][whpl]
        xid[i]['y'] = xid[i]['y'][whpl]
        xid[i]['sflux'] = xid[i]['sflux'][whpl]
        xid[i]['serr'] = xid[i]['serr']

        if savemap:
            outfile = config.CLUSSBOX + 'clus_get_xid_model_' + maps[i]['band'] + '.fit'
            writefits(outfile, data=model, header_dict=maps[i]['shead'])

    return xid
