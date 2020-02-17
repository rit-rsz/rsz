################################################################################
# NAME : test.py
# DATE STARTED : Feb 17, 2020
# AUTHORS : Benjamin Vaughan
# PURPOSE : This is a python version of Adi and Alfredo's IDL lensing software
# ported from DF_LensCat_FP.pro
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
from astropy.cosmology import LambdaCDM as LCDM

#these are some important constants
hb  = 70                     #Hubble constante
OMm = 0.3                    #Omega Matter (dark + baryonic)
OMl = 0.7                    #Omega Lambda
OMk = 1. - OMm - OMl         #Omega curvature
DH = 2997.92/hb              #Hubble Distance in Mpc => DH = (c * 1e-3)/Ho

def lenscat(map_size, c_z, pixsize=0.06, searchlensed=0.75, SLregime=(2,2)):
    '''
    Inputs:
    map_size : size of the input map
    c_z : the cluster's redshift, Lens_z in Alfredo's code
    pixsize : the size of each pixel in arcsec/pixel, pixscale in Alfredo's code
    searchlensed :
    SLregime : Area where we want to apply strong lensing, in arcmin
    '''

    #Parameters to do lensing outside the SL regime --------
    if SLregime[0] > map_size[0]:
        SLregime[0] = map_size[0]
    if SLregime[1] > map_size[1]:
        SLregime[1] = map_size[1]

    SL_pix = SL_regime * 60. / pixsize #converting strong lensing regime from arcmin to pixels

    #for the weakly lensed region we don't estimate the whole image plane but limit to a WLmask x WLmask area around the lensed source
    WLmask = 2.0 # reference scale to estimate WL region in arcmin
    WLmask = WLmask*60./pix_scale #WLmask now in pixels

    #IDL version for calculating the luminosity distance / (1 + len_z)**2 @  https://idlastro.gsfc.nasa.gov/ftp/pro/astro/lumdist.pro
    #DL1 = lumdist(Lens_z, H0 = hb, Omega_M = OMm, Lambda0 = OMl,/SILENT) / (1.0+Lens_z)^2, this is the call to lumdist in IDL

    #Python has a version in python @ https://docs.astropy.org/en/stable/cosmology/#functions
    #I assumed a LambdaCDM model based off of what was in Adi and Alfredo's code, but there are other options if that is the wrong model.

    #define a LambdaCDM model with an input hubble constant, matter content, and dark energy content
    cos_model = LCDM(H0=hb, Om0=OHm, Odec0=OMl)
    dl1 = cos_model.luminosity_distance(c_z) / (1 +c_z)**2

    #from what I understand in Alfredo's code c_z is the redshift of the cluster and I think main_source_z is the redshift of a center source ?
    #I'm not sure and there is code below written in IDL on what to do with the main_source_z parameter
    #applies a normalization to deflection.
    ###########################################################################################################################################
    '''
    ;if scaled to a certian source redshift (main_source_z):
    ;calculate angular diameter distances
    IF n_elements(main_source_z) EQ 1 THEN BEGIN

    	;estimate angular diameter distance to main_source_z
    	DS1 = lumdist(main_source_z, H0 = hb, Omega_M = OMm, Lambda0 = OMl,/SILENT) / (1.0+main_source_z)^2
    	;estimate distance between the lensing cluster and main_source_z
    	DLS1= (DS1*(1.+main_source_z) - DL1*(1.+Lens_z))/(1.+main_source_z)     ;WARNING: this is only valid for Omega_k = 0

    	;estimate the deflection field normalization
    	DLS1_DS1 = DLS1/DS1

    ENDIF ELSE DLS1_DS1 = 1. ;if main_source_z was not provided the deflection field normalization should (very likely) be 1
    '''
    ###########################################################################################################################################
    #for this next part we need ADi's READ_ZITRIN script, which I didn't see anywhere in the email
    '''
    ;we read the lensing deflection fields provided by Adi Zitrin, and which have npixels x npixels elements as a text file
    ;we use the function READ_ZITRIN which takes as input the number of pixels in one side of the map (npixels) and the name of the file
    	wl_alpha_x = READ_ZITRIN(npixels,FILE_DFX)     ;read the deflection field in X
    	wl_alpha_y = READ_ZITRIN(npixels,FILE_DFY)     ;read the deflection field in Y

    	wl_alpha_x = TRANSPOSE(wl_alpha_x)             ;need to transpose (at least in IDL)
    	wl_alpha_y = TRANSPOSE(wl_alpha_y)             ;need to transpose (at least in IDL)
    '''
    ###########################################################################################################################################
