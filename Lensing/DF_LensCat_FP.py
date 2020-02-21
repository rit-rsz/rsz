################################################################################
# NAME : test.py
# DATE STARTED : Feb 17, 2020
# AUTHORS : AMB, Benjamin Vaughan
# PURPOSE : This is a python version of Adi and Alfredo's IDL lensing software
# ported from DF_LensCat_FP.pro
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :  input_cat: input catalogue to be lensed.
#       FILE_DFX / FILE_DFY: deflection field in X & Y (provided by Adi Zitrin. Use READ_ZITRIN.pro read this files).
#       Lens_z: cluster (lens) redshift.
#       Pix_Scale: pixel size (Pix_Scale x Pix_Scale) in arcsec (default: 0.06" x 0.06").
#       LL: size of square map (L x L) in arcmin.
#       NPIXELS: number of pixels in map (NPIXELS).
#       SEARCHLENSED: this number is to find the lensed source which can be extended (arcs). Use the default value 0.75
#
# OUTPUTS : an array with [source ID (input), X position (input), Y position (input), redshift (input),
# 				unlensed flux (input), lensed X position, lensed Y position, magnification]
# 		Note that, to build the lensed map, you have to you should use the lensed positions and unlensed flux*magnification
#
# REVISION HISTORY :
################################################################################
from astropy.cosmology import LambdaCDM as LCDM
import meshgrid.py
import read_zitrin.py

#these are some important constants
# Should we use the origianl 0.7 for hb and do Ho = 100*hb to stay consistant with the original code
hb  = 0.7                    #Hubble constant (Ho = 100*hb)
OMm = 0.3                    #Omega Matter (dark + baryonic)
OMl = 0.7                    #Omega Lambda
OMk = 1. - OMm - OMl         #Omega curvature
# Check value now that we are using hb = 70 instead of 0.7
DH = 2997.92/hb              #Hubble Distance in Mpc => DH = (c * 1e-3)/Ho

'''
using Planck 2015 cosmology ... Ade+2015 a (Table 4 last column)
	hb  = 0.6774                    ; Hubble constante (Ho = 100*hb)
	OMm = 0.3089                    ; Omega Matter (dark + baryonic)
	OMl = 0.6911                    ; Omega Lambda
	OMk = 1.0 - OMm - OMl           ; Omega curvature
	DH  = 2997.92/hb                  ; Hubble Distance in Mpc => DH = (c * 1e-3)/Ho
'''

def lenscat(map_size, c_z, pixsize=0.06, searchlensed=0.75, SLregime=[2,2,2],DefField_Norm=None):
    '''
    Inputs:
    map_size : size of the input map
    c_z : the cluster's redshift, Lens_z in Alfredo's code
    pixsize : the size of each pixel in arcsec/pixel, pixscale in Alfredo's code
    searchlensed : default = 0.75 , find the lensed source to extend
    SLregime : Area where we want to apply strong lensing, in arcmin
    '''

    #Parameters to do lensing outside the SL regime --------
    ''' We will probably only have one map at a time, so don't need these three...'''
    if SLregime[0] > map_size[0]:
        SLregime[0] = map_size[0]
    if SLregime[1] > map_size[1]:
        SLregime[1] = map_size[1]
    if SLregime[2] > map_size[2]:
        SLregime[2] = map_size[2]

    SL_pix = SL_regime * 60. / pixsize #converting strong lensing regime from arcmin to pixels

    #for the weakly lensed region we don't estimate the whole image plane but limit to a WLmask x WLmask area around the lensed source
    WLmask = 2.0 # reference scale to estimate WL region in arcmin
    WLmask = WLmask*60./pix_scale #WLmask now in pixels

    #IDL version for calculating the luminosity distance / (1 + len_z)**2 @  https://idlastro.gsfc.nasa.gov/ftp/pro/astro/lumdist.pro
    #DL1 = lumdist(Lens_z, H0 = hb, Omega_M = OMm, Lambda0 = OMl,/SILENT) / (1.0+Lens_z)^2, this is the call to lumdist in IDL

    #Python has a version in python @ https://docs.astropy.org/en/stable/cosmology/#functions
    #I assumed a LambdaCDM model based off of what was in Adi and Alfredo's code, but there are other options if that is the wrong model.

    #define a LambdaCDM model with an input hubble constant, matter content, and dark energy content
    cos_model = LCDM(H0=hb*100.0, Om0=OMm, Odec0=OMl)
    '''
    Alternative cos model is directly using built in Planck15
    '''
    dl1 = cos_model.luminosity_distance(c_z) / (1 +c_z)**2

    #from what I understand in Alfredo's code c_z is the redshift of the cluster and I think main_source_z is the redshift of a center source ?
	#the luminosity distance functiom from astropy should be double checked
    #applies a normalization to deflection.

    if main_source_z == 1 and OMk == 0:
        #estimate angular diameter distance to main_source_z
        DS1 = cos_model.luminosity_distance(main_source_z) / (1.0+main_source_z)**2
        #estimate distance between the lensing cluster and main_source_z
        DLS1 = (DS1*(1. + main_source_z) - DL1*(1 + Lens_z)) / (1. + main_source_z) # WARNING: only true for Omega_k = 0

        #estimate the deflection field normalization
        DLS1_DS1 = DLS1 / DS1

    elif DefField_Norm != None: # DefField_Norm usually comes from another script.. here being passed as variable
        DLS1_DS1 = DefField_Norm

    else:
        DLS1_DS1 == 1. #if main_source_z was not provided the deflection field normalization should (very likely) be 1

    ###########################################################################################################################################

    '''
    # make sure the deflection fields are NxN
    if len(wl_alpha_x[:,0]) > len(wl_alpha_x[0,:]):
        wl_alpha_x = wl_alpha_x[0:len(wl_alpha_x[:,0])-1, :]
        wl_alpha_y = wl_alpha_x[0:len(wl_alpha_y[:,0])-1, :]

	elif len(wl_alpha_x[:,0]) < len(wl_alpha_x[0,:]):
		wl_alpha_x = wl_alpha_x[:, 0:len(wl_alpha_x[:,0])-1]
		wl_alpha_y = wl_alpha_x[:, 0:len(wl_alpha_x[:,0])-1]

	wl_lengthx = len(wl_alpha_x[:,0])
    wl_lengthy = len(wl_alpha_y[0,:])
    # Print check
    print('This should be 1 ...', wl_lengthx/wl_lengthy)

    YXarray = np.meshgrid(np.arange(wl_lengthy) + 1.0, np.arange(wl_lengthx) + 1.0)

    # Not sure which index is needed see below
    WLX = YXarray[:,:,0]
    WLY = YXarray[:,:,0]

    # Clear some memory
    YXarray = 0

    # this is an example of what meshgrid does
    # a =  MESHGRID(findgen(3),findgen(3))
    # print, a[*,*,0]
    #       0.00000      1.00000      2.00000
    #       0.00000      1.00000      2.00000
    #       0.00000      1.00000      2.00000
    # print, a[*,*,1]
    #       0.00000      0.00000      0.00000
    #       1.00000      1.00000      1.00000
    #       2.00000      2.00000      2.00000

    # Weak lensing: estimate gradients
    wl_dax_dy = np.arange(wl_lengthx,wl_lengthy)
    wl_dax_dx = np.arange(wl_lengthx,wl_lengthy)
    wl_day_dy = np.arange(wl_lengthx,wl_lengthy)
    wl_day_dx = np.arange(wl_lengthx,wl_lengthy)

    # Scipy has something that can maybe do this with misc derivative
    if wl_lengthx == wl_lengthy:
        for i in range(wl_lengthx):
            wl_dax_dy = diff(wl_alpha_x[i:])
            wl_dax_dx = diff(wl_alpha_x[:i])
            wl_day_dy = diff(wl_alpha_y[i:])
            wl_day_dx = diff (wl_alpha_y[:i])

    else:
        for i in range(wl_lengthx):
            wl_dax_dx[:i] = diff(wl_alpha_x[:i])
            wl_day_dx[:i] = diff(wl_alpha_y[:i])

        for i in range(wl_lengthy):
            wl_dax_dy[i:] = diff(wl_alpha_x[i:])
            wl_day_dy[i:] = diff(wl_alpha_y[i:])

##################################################################################
    '''
    IF wl_lengthx EQ wl_lengthy THEN BEGIN ;at this point wl_lengthx should be equal to wl_lengthy ... but you never know
    	FOR ii=0., wl_lengthx-1. DO BEGIN
    		wl_dax_dy[ii,*] = DERIV(wl_alpha_x[ii,*])       ;deflection field X with respect to y
    		wl_dax_dx[*,ii] = DERIV(wl_alpha_x[*,ii])       ;deflection field X with respect to x
    		wl_day_dy[ii,*] = DERIV(wl_alpha_y[ii,*])       ;deflection field Y with respect to y
    		wl_day_dx[*,ii] = DERIV(wl_alpha_y[*,ii])       ;deflection field Y with respect to x
    	ENDFOR
    ENDIF ELSE BEGIN
    	FOR ii=0., wl_lengthx-1. DO BEGIN
    		wl_dax_dx[*,ii] = DERIV(wl_alpha_x[*,ii])       ;deflection field X with respect to y
    		wl_day_dx[*,ii] = DERIV(wl_alpha_y[*,ii])       ;deflection field X with respect to x
    	ENDFOR
    	FOR ii=0., wl_lengthy-1. DO BEGIN
    		wl_dax_dy[ii,*] = DERIV(wl_alpha_x[ii,*])       ;deflection field Y with respect to y
    		wl_day_dy[ii,*] = DERIV(wl_alpha_y[ii,*])       ;deflection field Y with respect to x
    	ENDFOR
    ENDELSE
    '''
##################################################################################


# Strong lensing arrays -------------------------------------------------------------------
# here I make a copy of the different arrays (deflection fields, gradients and also the meshgrid arrays)
# but only for the region were we will perform strong lensing analysis. This made thing easier for me when
# I was writting the code and doing tests
    # This seems like it might not be the right python syntax
    sl_alpha_x = wl_alpha_x[(wl_lengthx-SLregime_pix)/2.:(wl_lengthx+SLregime_pix)/2.,(wl_lengthx-SLregime_pix)/2.:(wl_lengthx+SLregime_pix)/2.]
    sl_alpha_y = wl_alpha_y[(wl_lengthx-SLregime_pix)/2.:(wl_lengthx+SLregime_pix)/2.,(wl_lengthx-SLregime_pix)/2.:(wl_lengthx+SLregime_pix)/2.]

    sl_dax_dy = wl_dax_dy[(wl_lengthx-SLregime_pix)/2.:(wl_lengthx+SLregime_pix)/2.,(wl_lengthx-SLregime_pix)/2.:(wl_lengthx+SLregime_pix)/2.]
    sl_dax_dx = wl_dax_dx[(wl_lengthx-SLregime_pix)/2.:(wl_lengthx+SLregime_pix)/2.,(wl_lengthx-SLregime_pix)/2.:(wl_lengthx+SLregime_pix)/2.]
    sl_day_dy = wl_day_dy[(wl_lengthx-SLregime_pix)/2.:(wl_lengthx+SLregime_pix)/2.,(wl_lengthx-SLregime_pix)/2.:(wl_lengthx+SLregime_pix)/2.]
    sl_day_dx = wl_day_dx[(wl_lengthx-SLregime_pix)/2.:(wl_lengthx+SLregime_pix)/2.,(wl_lengthx-SLregime_pix)/2.:(wl_lengthx+SLregime_pix)/2.]

    sl_lengthx = len(sl_alpha_x[:0])

    YXarray = np.meshgrid(np.arrange(sl_lengthx)+1, np.arrange(sl_lengthx)+1)
    SLX = YXarrray[::0]
    SLY = YXarray[::1]
    YXarray = 0

    # shift between WL and SL frames (you will need this number to estimate the correct positions of the weak and strong lensed sources)
    pos_shift = (wl_lengthx-sl_lengthx)/2.

    #------------------------------------------------------------- OJO -----------------------------------------------------------------

    # Read source catalogue (input_cat)
    read_default_format = 'l,d,d,f,f'   #Default Format
    print('Reading Catalogue to be lensed')
    # Uses read call but this kind of strange. This is I think reading lines as varriables like dict keys
    # readcol, input_cat, dfid, dfxx, dfyy, dfredz, dfflux, skipline = 1, FORMAT = read_default_format

    # Convert coordinates from arcsecs to pixels in the larger WL frame
    xx = (dfxx + LL * 60./2.) / Pix_Scale
    yy = (dfyy + LL * 60./2.) / Pix_Scale
    # Number of sources in catalogue
    Nsources = len(xx)


    for i in range(NSources):
    # Some code here for timeing that we can skip
        '''
   ti = systime(/seconds)          ; don't worry about this, it's only to know when the lensing analysis started
   PRINT, 'Lensing source = ', ii+1, ' / ', NSources      ; don't worry about this, it just prints which source you are dealing with
   '''
# begin dealing with foreground sources (we don't need to lens these) ---------------------

	    if len_z > dfredz:
	        print('Source ',i+1,' is a foreground source')
	        mindX = xx[i]
	        mindY = yy[i]
	        mmu = 1.0
	    else:
	    # Finish dealing with foreground sources --------------------------------------------------

	    # Begin dealing with background sources (we do need to lens these)
	        # Estimate angular diameter distance to source (z = dfredz)
	        DS2 = cos_model.luminosity_distance(dfredz[i]) / (1.0+dfredz[i])**2
	        # Estimate distance between the lensing cluster and source
	        DLS2 = (DS2*(1.+dfredz[i]) - DL1*(1.+Lens_z)
			#this is the factor to scale the deflection field to different source redshifts.
			scaling_factor = (DLS2/DS2) / DLS1_DS1


	# ; begin SL regime analysis
				# print('time 1 = %s',                ;don't worry about this
	# ;			IF SQRT(dfxx[ii]^2.+dfyy[ii]^2.) LE SLregime*60. THEN BEGIN    ;don't worry about this

				#check if source is within strong lensing region, need to have a better understanding of what the 60 and 2 mean?
				#2 is since we are looking at the center of the strong lensing regime and that gives us the distance to an edge
				#60 is a conversion to arc minutes ?
				if abs(dfxx[i]) < SLregime[0] * 60. / 2. and abs(dfyy[i]) < SLregime[1] * 60 / 2.:
					print('Starting on SL regime')
					#estimate source position in SL frame
					source_x = (dfxx[i] + SLregime * 60 / 2.) / pixsize #we might need to add 1 here -Alfredo
					source_y = (dfyy[i] + SLregime * 60 / 2.) / pixsize #we might need to add 1 here -Alfredo
					#estimate the magnifications (from lensing eq.), adding sl_day_dx and sl_day_dy element by element
					poission = np.add(sl_dax_dx, sl_day_dy) * scaling_factor
					magnification = abs(1 / (1 - poission + (np.matmul(sl_dax_dx, sl_day_dy) - np.matmul(sl_day_dy,sl_dax_dx)) * scaling_factor))
					#will need to double check that these matrix operations work on sl_day_dx and sl_dax_dx
					#find the pixels where lensed image (or images) end up
					source_dist = sqrt( (SLX-sl_alpha_x*scaling_factor - source_x)**2 + (SLY-sl_alpha_y * scaling_factor - source_y)**2)
					indXY = np.where(source_dist < searchlensed) #this may have to change for a 2D array

					#if a lensed source was found within searched lensed do this
					if len(indXY) > 0:
						#cannot find an ij function in IDL anywhere but this is supposed to find indices for i and j corresponding to image pixels
						indXY = ij( indXY, sl_lengthx)

						#here we check for multiplicity
						#if more than 1 pixel satisfy the source_dst < searchlensed arg we might have multiplicty
						if len(indXY) > 1:
							#cut a square sub-map including all pixels with "source_dist < searchlensed" for easy computation
							min_x = min(indXY[:,0])
							min_y = min(indXY[:,1]) #not sure here i feel like this should be [0,:]
							max_x = max(indXY[:,0])
							max_y = max(indXY[:,1]) #ditto

							temp_multi_x = max_x - min_x
							temp_multi_y = max_y - min_y
							if temp_multi_x - temp_multi_y >= 0:
								multi_ll = temp_multi_x
							else:
								multi_ll = temp_multi_y

							#there are a bunch of if cases to consider here
							if min_x + multi_ll) < len(source_dist[:,0]):
								if min_y + multi_ll < len(source_dist[0,:]):

									regmap = source_dist[min_x:min_x + multi_ll, min_y:min_y+multi_ll]
								else:
									regmap = source_dist[min_x:min_x + multi_ll, min_y-1:len(source_dist[:,0])-1]
							elif min_y + multi < len(source_dist[0,:]):
								regmap = source_dist[min_x-1:len(source_dist[:,0]), min_y:min_y + multi_ll]
							else:
								regmap = source_dist[min_x-1:len(source_dist[:,0])-1, min_y-1:len(source_dist[0,:])-1]

							#indXY2 -the second one (we should come up with a better name maybe)
							indXY2 = np.where(regmap < searchlensed) #again, still unsure how this will pan out for 2D arrays

							indXY2 = ij( indXY2, len(regmap[:,0])) #still can't find documenation for an ij function in idl
							#don't understand what this is supposed to be doing aside from surface level finding i,j indices
							#corresponding to image pixels
							reg_centroids = np.zeros(2, len(indXY[:,0])) #create an empty array to put stuff in

							for j in range(len(indXY2[:,0])):
								#----------------- some idl code here
								#region = SEARCH2D(regmap, indXY2[jj,0], inXY2[jj,1], 0, searchlensed, /DIAGONAL)
								#regmap is the input map
								#indXY2[j,0] is initial x position
								#indXY2[j,1] is initial y position
								#0 is the lower limit
								#searchlensed is the upper limit
								#link to documentation for SEARCH2D https://www.harrisgeospatial.com/docs/SEARCH2D.html

								#initial an empty array structure to put stuff in!
								regmask = np.zeros(len(regmap[:,0]), len(regmap[:,0]))
								#give a value of 1 to all pixels in our region !
								regmask[region] = 1.0

								#find the center of mass with mask applied !
								reg_centroids[:, j] =
