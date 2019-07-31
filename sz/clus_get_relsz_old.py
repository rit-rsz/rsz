
################################################################################
# NAME : clus_get_relsz.py
# DATE STARTED : June 18, 2019
# AUTHORS : Dale Mercado
# PURPOSE : Given some nominal SZ inputs, this returns the expected SZ
#      effect signal.  Output is in W m^-2 Hz^-1 sr^-1.
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#        nu  - the range of frequncies requested [fltarr, GHz]
#        y   - the Compton y parameter, usually O(10^-4) [float, unitless]
#        T_e - the electron temperature of the gas [float, keV
#
#
# OUTPUTS :
# NOTES :
#       Code is based on calculations in these papers:
#       http://adsabs.harvard.edu/abs/1998ApJ...499....1C
#       http://adsabs.harvard.edu/abs/1998ApJ...502....7I
#       http://adsabs.harvard.edu/abs/2000ApJ...536...31N
#       http://adsabs.harvard.edu/abs/2002ARA%26A..40..643C
# REVISION HISTORY :
################################################################################

# Currently this is in a rough state, while I am trying to interperate the idl script
import scipy.io
import numpy as np
# from config import * #(this line will give me access to all directory variables)
from math import *
from astropy.io import fits
from scipy.interpolate import interp1d


def get_aij():

    aij = [[ 4.212483e-3, -4.444816e-1, 1.073238e+1, -1.096539e+2, 5.964685e+2, -1.925063e+3, 3.876921e+3, -4.927770e+3, 3.843121e+3, -1.677474e+3, 3.131417e+2],\
    [ -6.65860e-2, 3.141377e+0, -5.060928e+1, 2.931651e+2, -4.799553e+2, -1.656944e+3, 9.022748e+3, -1.777832e+4, 1.807900e+4, -9.401495e+3, 1.968976e+3],\
    [ 1.064073e+0, -1.495849e+1, 1.862757e+2, -1.026793+3, 2.276656e+3, -1.759309e+3, 1.512946e+3, -8.627984e+3, 1.693099e+4, -1.283544e+4, 3.354941e+3],\
    [ -1.110801e+1, 5.585558e+1, -3.339812e+2, 1.705455e+3, -2.884721e+3, -8.599791e+2, 7.487214e+3, -7.967433e+3, -3.628501e+2, 7.331197e+3, -4.181022e+3],\
    [ 6.567949e+1, -2.018768e+2, 1.693840e+2, -1.350618e+3, 2.741178e+3, -1.865475e+3, 7.054196e+3, -1.223620e+4, -4.349851e+3, 2.464867e+4, -1.444005e+4],\
    [ -2.285280e+2, 5.817376e+2, 6.539315e+2, -9.223397e+2, 3.300752e+2, -3.782007e+3, 5.413227e+3, -8.423819e+3, -4.543510e+3, 2.241957e+4, -1.262084e+4],\
    [ 4.882355e+2, -1.136911e+3, -1.737435e+3, 2.902696e+3, -3.374909e+2, 1.583587e+3, 4.285821e+3, -2.088839e+3, -8.443411e+3, 1.643675e+4, -7.819763e+3],\
    [ -6.493386e+2, 1.434664e+3, 1.998541e+3, -2.732035e+3, 2.384819e+1, -5.098857e+3, 2.130640e+3, 1.173659e+3, -6.531589e+3, 2.654111e+3, -1.648207e+3],\
    [ 5.247520e+2, -1.134319e+3, -1.008644e+3, 1.167943e+2, 3.264981e+3, -2.006417e+3, 3.974758e+3, 4.929594+3, -7.369537+3, 1.842694e+3, 3.854699e+3],\
    [ -2.360646e+2, 5.147089e+2, 1.003613e+2, 6.701438e+2, -1.143976e+3, 1.427423e+3, -8.377902e+3, 7.569877e+3, -2.057870e+2, -5.147649e+3, 1.635543e+3],\
    [ 4.537578e+1, -1.024535e+2, 4.040647e+1, 1.039578e+2, -2.327177e+2, 7.547122e+3, -1.097034e+4, 9.269639e+3, -3.302366e+3, 4.663621e+2, -2.266015e+2]]

    return aij


def get_Y(x):

    xp = x * coth(x/2)
    sp = x / sinh(x/2)

    y0 = -4 + xp

    y1 = -10 + (47/2) * xp - \
        (42/5) * xp**2 + (7/10) * xp**3 + \
        sp**2 * ((-21/5) + (7/5)*xp)

    y2 = (-15/2) + (1023/8) * xp - (868/5) * xp**2 + (329/5) * xp**3 - \
        (44/5) * xp**4 + (11/30) * xp**5 + \
        sp**2 * (-434/5 + (658/5) * xp - \
        (242/5) * xp**2 + (143/30) * xp**3) + \
        sp**4*(-44/5 + (187/60) * xp)

    y3 = (15/2) + (2505/8) * xp - (7098/5) * xp**2 + (14253/10) * xp**3 - \
        (18594/35) * xp**4 + (12059/140) * xp**5 - (128/21) * xp**6 + \
        (16.0d0/105) * xp**7 + \
        sp**2 * ((-7098/10) + (14253/5) * xp - (102267/35) * xp**2 + \
        (156767/140) * xp**3 - (1216/7) * xp**4 + (64/7) * xp**5) + \
        sp**4*((-18594/35) + (205003/280) * xp - (1920/7) * xp**2 + \
        (1024/35) * xp**3) + sp**6 * ((-544/21) + (992/105) * xp)

    y4 = (-135/32) + (30375/128) * xp - (62391/10) * xp**2 + \
        (614727/40) * xp**3 - (124389/10) * xp**4 + (355703/80) * xp**5 - \
        (16568/21) * xp**6 + (7516/105) * xp**7 - (22/7) * xp**8 + \
        (11/210) * xp**9 + \
    sp**2 * ((-62391/20) + (614727/20) * xp - (1368279/20) * xp**2 + \
    (4624139/80) * xp**3 - (157396/7) * xp**4 + \
    (30064/7) * xp**5 - (2717/7) * xp**6 + (2761/210) * xp**7) + \
    sp**4 * ((-124389/10) + (6046951/160) * xp - (248520/7) * xp**2 + \
    (481024/35) * xp**3 - (15972/7) * xp**4 + (18689/140) * xp**5) + \
    sp**6 * ((-70414/21) + (465992/105) * xp - (11792/7) * xp**2 + \
    (19778/105) * xp**3) + sp**8*((-682/7) + (7601/210) * xp)

    # This should be a dictionary
    bigY = {'y0':y0,'y1':y1,'y2':y2,'y3':y3,'y4':y4}

    return bigY




def clus_get_relsz(nu,y,T_e,\
                        canonical = 0,x = 0,nor = 0,\
                        szpack=1,recomplookup=0,\
                        lookup=1,dIfidp=0,header = 0)
                        # in idl ones that didnt have an input were
                          # canonical,x,nor,difid,header

    errmsg = False
    # If there is any dIfidp input set it equal to 1
    if dIfidp:
        dIfidp = 1

    T_0     = 2.725                   #CMB temperature, K
    k_B     = 1.3806503e-23           #Boltzmann constant, J/K
    h       = 6.626068e-34            #Planck constant, J s
    m_e     = 5.11e2                  #electron mass keV
    c       = 2.99792458e8            #m/s
    HztoGHz = 1e9                     #Hz -> GHz
    mks_to_MJysr = 1e20               #get me to MJy/sr!

    I0 = 2 * (k_B * T_0)**3 / (h * c)**2

    if recomplookup:
        #this does some sxaddpar header stuff but doesnt return anythin and has a stop?
        #Currently is not written recomlookup set zero default so skipping for now
        clus_szdi_lookup()

    if lookup:

        deltaI = np.zeros(nu.size)
        if dIfidp = 0:
            dIfid = fits.open('lookup/clus_szdi_lookup.fits')

        npts = dIfit[0].header['npts']
        dy = dIfit[0].header['dy']
        dT  =dIfit[0].header['dT']

        # Is this still necessary??
        # this unnecessary replicate is to fool interpolate, which is dumb
        thisy = npts * y / 1e-3
        thisx = npts * T_e / 25

        if thisy > npts:
            errmsg = 'Input y is larger than test grid.'
            return None, errmsg
            # Talk with Ben about returning this ERRMSG
        if thist > npts:
            errmsg = 'Input T_e is larger than test grid.'
            return None, errmsg

        for i in range(nu):
            if abs(nu[i] - 140.187) < 0.1:
                knu = 0
            if abs(nu[i]) - 600) < 0.1:
                knu = 1
            if abs(nu[i]) - 857.143) < 0.1:
                knu = 2

            if knu.size = 0:
                errmsg = 'Confused about nu mapping'
                return None, errmsg
                # Again talk with Ben about the ERRMSG
            # This gets reformed and then interpolated??? use interpld or use Rbf in scipy
            # Need to better understand reform
            # deltaI[i] = interp1d(reshape(dIfid[:,:,knu]),thisy,thist)

        # This 1 is apparently
        print('1')
        # Need to fix something with the syntax of this retun. the else only activates if there is a : afer the return???
        return deltaI,errmsg

    else:

#       szpack doesnt work at 0 electron temperature
        if abs(T_e) < 1:
                                # T_e = 0
                                # canonical = 1
        szpack = 0

        if szpack:
            tau = y * 420 / T_e**0.9
            te = T_e

            if not te:
                te = 10.0
            else:
                te = te

            if not vpec:
                vpec = 0
            else:
                vpec = vpec

            if not ngrid:
                ngrid = 100
            else:
                ngrid = ngrid

            # it now opens and writes to lookup/szparams.txt

            # Lines 237-277 are stange I need to make sure that I know exactly what they are doing as well as that I am doing them right

        else:

#           check that the nu input isn't crazy
            if min(nu) < 1e-4 or max(nu) > 1e4:
                errmsg = 'Frequency input out of bounds, please check : ' + errmsg
                return None, errmsg

            y = np.double(y)
#           Check to make sure the y value isn't crazy
            if y < 1e-7 or y > 1e-2:
                errmsg = 'T_e input out bounds, please check : ' + errmsg
                return None, ermmsg

            T_e = np.double(T_e)
#           Check that the T_e input isn't beyond the equation
            if T_e < 0 or T_e > 200:
                errmsg = 'T_e input out of bounds, please check : ' + errmsg
                return None, errmsg

            if T_e = 0:
                canonical = 1

#           override T_e if canonical is set
            if canonical:
                T_e = 0
                if not nor:
                    nor = 0

#           calibrate nu up to something easier to use
            nup = nu * HztoGHz

#           make some of the fundamental frequency arrays
            nx = len(nu)
#           note here that X can be an output!
            X = (h * nup) / (k_b * T_0)
            Z = (1 / 17.6) * (X - 1 * 2.4)

#           make the scaled electron temperatures
            theta_e = T_e / m_e

#           here we're sticking hard to Nozawa's upper paramter limits;
#           there will be an error if they're out of bounds
            if theta_e > 2e-2 and not nor:
                if theta_e < 5.1e-1:
                    Theta = 25 * (theta_e - 0.01)
#                   get the a_ij we need for the corrections
                    aij = get_aij()
#                   and make an ancillary stuff we need for the matric operation
                    nterms = aij.size
                    power = np.zeros(0,nx)
#                   loop through all frequencies
#                   IDL marks subtracted from 1L but I need to check if this is they way to do it in python
                    for ix in range(nx-1L):
                        if X[ix] > 2.5:
                            Z_vec = Z[ix]**(np.transpose(power))
                            R[ix] = total(aij * Theta_vec * Z_vec)
#                   if the T_e is too high (should have caught that already)
#                   bail with an error
                else:
                    errmsg = \
                        'T_e caused calculation to go out of bounds, please check : ' + \
                        errmsg
                    return None, errmsg

            else:
                R = 0

            if T_e = 0:
                theta_e = 0

#           get the big Y terms
            bigYs = get_Y(X)
#           compute there sum; note this goes to the canonical SZ spectrum in
#           the limit T_e = 0
#           Need to check the syntax on bigYs.y0
#           Need to look up this math to check order of oporations
            Yterm = bigYs.y0 + theta_e * bigYs.y1 + theta_e**2 * bigYs.y2 + \
                    theta_e**3 * bigYs.y3 + theta_e**4 * bigYs.y4
#           Make Nozawa's F(theta_e,X)
            FoftX = Yterm * (X * exp(X) / (exp(X)-1)) + R * theta_e
            FoftY = bigYs.y0 * (X * exp(X) / (exp(X)-1))

            npts = X.shape

#           make dn / n_0 term
            dnn = FoftX * repeat(y,npts)

#           make the term involving X alone
            xterm = X**3 / (exp(X) - 1)

#           make I_0, the term that calibrates this all to W m^-2 Hz^-1
            I0 = repeat(2 * (k_B * T_0)**3 / (h*c)**2), npts)

#           and multiply them all together to get something useful
            deltaI = xterm * dnn * I0

#           if we got this far we propbably have a correct SZ effect shape in
#           our band and can get out of here
            return deltaI, errmsg

#   Error handler has been replaced by the returns

    return deltaI, errmsg
