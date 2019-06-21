################################################################################
# NAME : get_xid.py
# DATE STARTED : June 21, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE : idk lol
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
    cols = ['PSW', 'PMW', 'PLW']
    ncols = len(cols)
    psfs = np.empty(ncols)

    for i in range(ncols):
        if cols[i] in maps[i]['band']:
            index = i
        if i == 0:
            pswfits = maps[index]['file']
        elif i == 1:
            pmwfits = maps[index]['file']
        elif i == 2:
            plwfits = maps[index]['file']
        psfs[i] = get_spire_beam_fwhm(cols[i])

    xra = cat['ra']
    xdec = cat['dec']
    f_src = cat['flux']

    if len(xra) >= 2048:
        xra = xra[0:2048]
        xdec = xdec[0:2048]
        f_src = xdec[0:2048]
  # I do not know wwhat the following code does in IDL so its hard for me to write a comparison in python
  # IF N_ELEMENTS(xra) GT 2048 THEN BEGIN
  #    whpl = SORT(f_src)
  #    xra = xra[whpl[0:2048]]
  #    xdec = xdec[whpl[0:2048]]
  #    f_src = f_src[whpl[0:2048]]
  # ENDIF

    name = str(maps[0]['name']) + '_xid'

    astromcorr = [0.0, 0.0]

    if verbose:
        print('Calling XID')

    # don't know what the below function does at all, can't find anything on internet.
      # LSTDRV_READSMAP,pswfits,pmwfits,plwfits,xra,xdec,psfs,spirecat,$
      #                 f_src=f_src,name=name,fcatmod=spirecat2,$
      #                 astromcorr=astromcorr,thresh=thresh,noconf=~confusionerrors,noseg=1

    if verbose:
        print('XID done')

    xid = []

    for i in range(ncols):
        #where did spire cat come from
        sid = spirecat2[i]['xid']
        sra = spirecat2[i]['inra']
        sdec= spirecat2[i]['indec']
        pflux = spirecat2[i]['obs250'] / mJy2Jy
        perr = spirecat2[i]['e250'] / mJy2Jy
        sflux = spirecat2[i]['obs250'] / mJy2Jy
        serr = spirecat2[i]['e250'] / mJy2Jy

    #call to convert from RA/DEC to other coordinates...

    #can't find an equivalent for this.
   # ;; make a fake image from the xid catalog
   # model = IMAGE_MODEL(x,y,sflux,(*maps[icol]).astr.naxis[0],$
   #                   (*maps[icol]).astr.naxis[1],(*maps[icol]).psf)

        whpl = np.where(spirecat[i]['obs250'] / mJy2Jy >= 0.0)

        xid_sb = {'id' : sid[whpl], 'band' : cols[i],
                  'ra' : sra[whpl], 'dec' : sdec[whpl],
                  'x' : x[whpl], 'y' : y[whpl],
                  'flux' : sflux[whpl], 'err' : serr[whpl],
                  'fluxi' : pflux[whpl], 'erri' : perr,
                  'model' : model, 'mf' : mf}

        xid.append(xid_sb) #i've had issues with code like this before not sure if this will work out.
        if savemap:
            outfile = config.CLUSBOS + 'clus_get_xid_model_' + cols[i] + '.fits'
            #writefits call !

    return xid, err

def sort(arr):
    n = len(arr)
    for i in range(n):
        for j in range(0, n-i-1):
            if arr[j] > arr[j+1]:
                temp = arr[j+1]
                arr[j+1] = arr[j]
                arr[j] = temp
    return arr
