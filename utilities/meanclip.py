################################################################################
# NAME : subtract_xcomps.py
# DATE STARTED : June 21, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE : computes an iterativley sigma-clipped mean on a data set.
# https://gist.github.com/jiffyclub/1310947/729dd001644a534b998503a9d8a81a0834cf8908
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################import numpy


def meanclip(indata, clipsig=3.0, maxiter=5, converge_num=0.02, verbose=0):
    skpix = indata.reshape( indata.size, )

    ct = indata.size
    iter = 0; c1 = 1.0 ; c2 = 0.0

    while (c1 >= c2) and (iter < maxiter):
        lastct = ct
        medval = numpy.median(skpix)
        sig = numpy.std(skpix)
        wsm = numpy.where( abs(skpix-medval) < clipsig*sig )
        ct = len(wsm[0])
        if ct > 0:
            skpix = skpix[wsm]

        c1 = abs(ct - lastct)
        c2 = converge_num * lastct
        iter += 1
    # End of while loop

    mean  = numpy.mean( skpix )
    sigma = robust_sigma( skpix )

    if verbose:
    prf = 'MEANCLIP:'
    print '%s %.1f-sigma clipped mean' % (prf, clipsig)
    print '%s Mean computed in %i iterations' % (prf, iter)
    print '%s Mean = %.6f, sigma = %.6f' % (prf, mean, sigma)

    return mean, sigma
