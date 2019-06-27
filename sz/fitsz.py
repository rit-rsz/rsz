################################################################################
# NAME : fitsz.py
# DATE STARTED : June 24, 2019
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

def fitsz(radave, params, beam=None, maxlim=3600, minlim=0, noweight=1, superplot=0, verbose=1):
    ncols = len(radave)

    increment = np.empty(ncols)
    offsets = np.empty(ncols)
    fit = np.empty(2, ncols)

    for i in range(ncols):
        rc = np.tile(params['rc'], radave[i]['midbin'].shape)
        beta = np.tile(params['beta'], radave[i]['midbin'].shape)
        if beam:
            for j in range(radave[i]['midbin'].shape[0]):
                for k in range(radave[i]['midbin'].shape[1]):
                    roft = (1.0 + (radave[i]['midbin'] / rc)**2)**((1.0 - 3.0 * beta) / 2.0)
        else:
            xrad = []
            for j in range(np.amax.(radave[i]['midbin']), 0, -1):
                xrad.append(j)
            xrad = np.array(xrad)
            xrc = np.tile(params['rc'],2*np.amax(radave[i]['midbin']))
            xbeta = np.tile(params['beta'], 2*np.amax(radave[i]['midbin']))
            roftprime = np.empty(len(radave[i]['midbin']))
            psf = np.empty(len(xrad))
            for j in range(len(radave[i]['midbin'])):
                roftprime[j] = (1.0 + (xrad[j] / xrc[j])**2)**((1.0 - 3.0 * xbeta[j]) / 2.0)
                pf = -1*(2.18 * math.log(2.0) / beam[i]**2)
                psf[j] = math.exp(pf * xrad[j]**2)
                psf[j] = psf[j] / np.sum(psf)
                #roftprimep = CONVOL(roftprime, psf, EDGE_TRUNCATE)
                roftprimep[j] = roftprimep[j] / np.amax(abs(roftprimep))
                #roft = INTERPL(roftprimep, xrad, radwave[i]['midbin'])

        whpl = []
        for j in range(radave[i]['fluxbin'].shape[0]):
            if np.isfinite(radave[i]['fluxbin'][j] == True) and radave[i]['midbin'][j] < maxlim:
                whpl.append(j)

        if superplot:
            #code to plot the stuffs this could probably be done in matplotlib...

        if noweight:
            #fitSZ = LINFIT(roft[whpl], radave[i]['fluxbin'][whpl], CHISQ, sigma)
        else:
            #fitSZ = LINFIT(roft[whpl], radave[i]['fluxbin'][whpl], MEASURE_ERRORS=radave[i]['errbin'][whpl], chisq=chisq, sigma=SIGMA)

        print(fitSZ)
        print(chisq) #probably gets made in LINFIT.
        print(sigma) #probably gets made in LINFIT.

        if superplot:
            #plot the new stuffs, again could probably be done in matplotlib, but i need to do real testing to find out.

        increment[i] = fitSZ[i]
        offsets[i] = fitSZ[0]
        fit[:, i] =fitSZ

    return fit
