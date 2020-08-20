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
#                 unlensed flux (input), lensed X position, lensed Y position, magnification]
#         Note that, to build the lensed map, you have to you should use the lensed positions and unlensed flux*magnification
#
# REVISION HISTORY :
################################################################################
from astropy.cosmology import LambdaCDM as LCDM
# from read_zitrin import *
import numpy as np
import sys
import matplotlib.pyplot as plt
from math import *
from astropy.wcs import WCS as world
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
from scipy.ndimage.measurements import label
from astropy.convolution import Gaussian2DKernel
import time
from lensing_utils import *


#these are some important constants
hb  = 0.7                    #Hubble constant (Ho = 100*hb)
OMm = 0.3                    #Omega Matter (dark + baryonic)
OMl = 0.7                    #Omega Lambda
OMk = 1. - OMm - OMl         #Omega curvature
DH = 2997.92/hb              #Hubble Distance in Mpc => DH = (c * 1e-3)/Ho

# using Planck 2015 cosmology ... Ade+2015 a (Table 4 last column)
#     hb  = 0.6774                    ; Hubble constante (Ho = 100*hb)
#     OMm = 0.3089                    ; Omega Matter (dark + baryonic)
#     OMl = 0.6911                    ; Omega Lambda
#     OMk = 1.0 - OMm - OMl           ; Omega curvature
#     DH  = 2997.92/hb                  ; Hubble Distance in Mpc => DH = (c * 1e-3)/Ho




def lense(catalogue, lens_z, searchlensed, file_dfx, file_dfy, LL, pixsize=0.25, SLregime=2, main_source_z=None, DefField_Norm=None):
	'''
	This is the lensing code developed by Adi and Alfredo and ported over to python. The purpose of this is to apply lensing to a
	simulation of point sources given some deflection shields.

	Inputs:
	catalogue (npy filename)- the catalogue of point sources to be lensed, should have locations in arcseconds and centered at (0,0)
	#format should be a rec array with keys of RA (arcs), DEC (arcs), Flux (Jy), Redshift.
	lens_z (float_)- the red shift of the cluster.
	searchlensed (float)- parameter to determine if lensed sources can be extended
	file_dfx (txt filename)- the deflection field in x
	file_dfy (txt filename)- the deflection field in y
	LL (float)- the area of the map in arcminutes (assumes a square)
	pixsize (float) - the pixelsize (arcseconds) of the deflection field
	SLregime (float) - the area of the strong lensing regime in arcminutes (assumed to be at center of map)
	main_source_z (float)- redshift used to create the deflection fields
	DefField_Norm (float)- a parameter used if main_source_z is not given (can be calculated with main_source_z if not provided)
	#if neither main_source_z or DefField_Norm are provided then DefField_Norm is assumed to be 1.
	'''

	# Define input cluster redshift, source redshift to which deflection field is scaled to, and the pixel
	# scale (both are input provided by Adi Zitrin). Often the deflection field will be scaled to DLS/DS=1.

	#if SLregime is bigger than the mapsize make the SLregime the mapsize
	if SLregime > LL:
	    SLregime = LL

	print('SLregime :',SLregime)
	SL_pix = SLregime * 60 / pixsize #converting strong lensing regime from arcmin to pixels
	#for the weakly lensed region we don't estimate the whole image plane but limit to a WLmask x WLmask area around the lensed source
	WLmask = 2.0 # reference scale to estimate WL region in arcmin
	WLmask = WLmask / pixsize * 60.0 #WLmask now in pixels

	print('WL region (in pixels):', WLmask)


	#estimate the angular distance to the lensing cluster (z = lens_z)
	cos_model = LCDM(hb, OMm, OMl)
	DL1 = cos_model.luminosity_distance(lens_z).value / (1 + lens_z)**2  # * 100 to match units in Adi's code.

	print('DL1 (luminosity distance) :',DL1)

	#If this is scaled to a certain source redshift (main_sourcez) then calculate angular diameter distances
	if main_source_z != None:
	    #estimate angular diameter distance to main_source_z
	    DS1 = cos_model.luminosity_distance(main_source_z).value / (1 + main_source_z)**2
	    #estimate distance between the lensing cluster and main_source_z
	    DLS1 = (DS1*(1. + main_source_z) - DL1*(1 + lens_z)) / (1. + main_source_z) # WARNING: only true for Omega_k = 0

	    #estimate the deflection field normalization
	    DLS1_DS1 = DLS1 / DS1

	elif DefField_Norm != None: # if main_source_z not provided but deflection field norm was given.
	    DLS1_DS1 = DefField_Norm

	else:
	    DLS1_DS1 = 1. #if main_source_z was not provided the deflection field normalization should (very likely) be 1

	#read in the deflection fields provided by Adi Zitrin, (npixel X npixel)
	wl_alpha_x = np.loadtxt(file_dfx) #deflection field of x
	wl_alpha_y = np.loadtxt(file_dfy) #deflection field of y

	# make sure the deflection fields are NxN
	if wl_alpha_x.shape[0] > wl_alpha_x.shape[1]:
	    wl_alpha_x = wl_alpha_x[0:wl_alpha_x.shape[1], :]
	    wl_alpha_y = wl_alpha_y[0:wl_alpha_X.shape[1], :]

	elif wl_alpha_x.shape[1] < wl_alpha_x.shape[0]:
	    wl_alpha_x = wl_alpha_x[:, 0:wl_alpha_x.shape[0]]
	    wl_alpha_y = wl_alpha_y[:, 0:wl_alpha_x.shape[0]]

	#weak lensing estimates
	wl_lengthx = len(wl_alpha_x[:,0])
	wl_lengthy = len(wl_alpha_y[0,:])
	# Print check
	print('This should be 1 ...', wl_lengthx/wl_lengthy)

	wl_x = np.arange(start=1, stop=wl_lengthx + 1)
	wl_y = np.arange(start=1, stop=wl_lengthy + 1)

	WLY, WLX = np.meshgrid(wl_x, wl_y)

	#weak lensing gradients
	wl_dax_dy = np.zeros((wl_lengthx, wl_lengthy))
	wl_dax_dx = np.zeros((wl_lengthx, wl_lengthy))
	wl_day_dy = np.zeros((wl_lengthx, wl_lengthy))
	wl_day_dx = np.zeros((wl_lengthx, wl_lengthy))


	if wl_lengthx == wl_lengthy:
	    for i in range(wl_lengthx):
	        wl_dax_dy[i,:] = deriv(wl_alpha_x[i,:])
	        wl_dax_dx[:,i] = deriv(wl_alpha_x[:,i])
	        wl_day_dy[i,:] = deriv(wl_alpha_y[i,:])
	        wl_day_dx[:,i] = deriv(wl_alpha_y[:,i])

	else:
	    for i in range(wl_lengthx):
	        wl_dax_dx[:,i] = deriv(wl_alpha_x[:,i])
	        wl_day_dx[:,i] = deriv(wl_alpha_y[:,i])

	    for i in range(wl_lengthy):
	        wl_dax_dy[i,:] = deriv(wl_alpha_x[i,:])
	        wl_day_dy[i,:] = deriv(wl_alpha_y[i,:])


	#Strong Lensing
	#make a copy of the grids but only over the SL region

	dim_minx = int((wl_lengthx-SL_pix)/2)
	dim_plusx = int((wl_lengthx+SL_pix)/2) + 1 #numpy end interval is exclusive.
	dim_miny = int((wl_lengthy-SL_pix)/2)
	dim_plusy = int((wl_lengthy+SL_pix)/2) + 1

	sl_alpha_x = wl_alpha_x[dim_minx:dim_plusx,dim_minx:dim_plusx]
	sl_alpha_y = wl_alpha_y[dim_miny:dim_plusy,dim_miny:dim_plusy]

	sl_dax_dy = wl_dax_dy[dim_miny:dim_plusy,dim_miny:dim_plusy]
	sl_dax_dx = wl_dax_dx[dim_minx:dim_plusx,dim_minx:dim_plusx]
	sl_day_dy = wl_day_dy[dim_miny:dim_plusy,dim_miny:dim_plusy]
	sl_day_dx = wl_day_dx[dim_minx:dim_plusx,dim_minx:dim_plusx]
	sl_lengthx = len(sl_alpha_x[:,0])
	sl_lengthy = len(sl_alpha_y[0,:])

	sl_x = np.arange(start=1, stop=int(sl_lengthx+1))
	sl_y = np.arange(start=1, stop=int(sl_lengthy+1))

	SLY, SLX = np.meshgrid(sl_x, sl_y)
	# shift between WL and SL frames (you will need this number to estimate the correct positions of the weak and strong lensed sources)
	pos_shift = (wl_lengthx-sl_lengthx)/2.

	print('pos shift between sl and wl [pixels]:',pos_shift)

	cat = np.load(catalogue, allow_pickle=True)
	print('Reading Catalogue to be lensed')
	dfxx = cat.item().get('RA')
	dfyy = cat.item().get('DEC')
	dfflux = cat.item().get('Flux')
	dfredz = cat.item().get('Redshift')

	# Convert coordinates from degrees to pixels in the larger WL frame
	xx = (dfxx + LL * 30)/ pixsize
	yy = (dfyy + LL * 30)/ pixsize

	# Number of sources in catalogue
	nsources = len(xx)
	print('number of sources to be lensed:', nsources)

	print('Set Up Complete beggining to lense sources')

	start = time.time()

	for i in range(nsources):
	    if not i % 1000:
	        print('on source %s / %s' % (str(i), str(nsources)))
	#     print('lensing source %s' % str(i + 1))
	    # begin dealing with foreground sources (we don't need to lens these) ---------------------
	    if lens_z > dfredz[i]:
	#         print('source %s is a foreground source' % str(i + 1))
	        mindX = xx[i]
	        mindY = yy[i]
	        mmu = 1.0

	    else:
	        # Finish dealing with foreground sources --------------------------------------------------
	        # Begin dealing with background sources (we do need to lens these)
	        # Estimate angular diameter distance to source (z = dfredz)
	        DS2 = cos_model.luminosity_distance(dfredz[i]).value / (1.0+dfredz[i])**2  # * 100 to match units in Adi's code
	        # Estimate distance between the lensing cluster and source
	        DLS2 = (DS2*(1.+dfredz[i]) - DL1*(1.+lens_z)) / (1 + dfredz[i])
	        #this is the factor to scale the deflection field to different source redshifts.
	        scaling_factor = (DLS2/DS2) / DLS1_DS1

	        #begin SL regime analysis
	        x_check = abs(dfxx[i]) < SLregime * 30
	        y_check = abs(dfyy[i]) < SLregime * 30
	        xy_check = x_check and y_check
	        if xy_check:
	#             print('beggining strong lensing regime')
	            #estimate source position in SL frame
	            source_x = int((dfxx[i] + SLregime * 30) / pixsize)
	            source_y = int((dfyy[i] + SLregime * 30) / pixsize)
	            #estimate the magnifications (from lensing eq.), adding sl_day_dx and sl_day_dy element by element
	            sl_dax_day_dxdy = np.multiply(sl_dax_dx, sl_day_dy)
	            sl_day_dax_dydx = np.multiply(sl_dax_dy, sl_day_dx)

	            poisson = np.add(sl_dax_dx, sl_day_dy) * scaling_factor
	            magnification = abs(1 / (1 - poisson + (sl_dax_day_dxdy - sl_day_dax_dydx) * scaling_factor))
	            #will need to double check that these matrix operations work on sl_day_dx and sl_dax_dx
	            #find the pixels where lensed image (or images) end up
	            x_side = SLX - sl_alpha_x * scaling_factor - source_x
	            y_side = SLY - sl_alpha_y * scaling_factor - source_y
	            source_dist = np.sqrt(x_side**2 + y_side**2)
	            indX, indY = np.where(source_dist < searchlensed) #this may have to change for a 2D array
	            #if a lensed source was found within searched lensed do this

	            if len(indX) > 1:

	                #cut a square sub-map including all pixels with "source_dist < searchlensed" for easy computation
	                min_x = min(indX)
	                min_y = min(indY)
	                max_x = max(indX)
	                max_y = max(indY)
	                temp_multi_x = max_x - min_x
	                temp_multi_y = max_y - min_y


	                #pick the larger of the two sides
	                if temp_multi_x - temp_multi_y >= 0:
	                    multi_ll = temp_multi_x
	                else:
	                    multi_ll = temp_multi_y



	                if (min_x + multi_ll) < len(source_dist[:,0]):
	                    if (min_y + multi_ll) < len(source_dist[0,:]):
	                        regmap = source_dist[min_x:min_x + multi_ll + 1, min_y:min_y+multi_ll + 1]
	                    else:
	                        regmap = source_dist[min_x:min_x + multi_ll+ 1, min_y:source_dist.shape[0] + 1]
	                elif (min_y + multi_ll) < len(source_dist[0,:]):
	                    regmap = source_dist[min_x:source_dist.shape[0]+1, min_y:min_y + multi_ll + 1]
	                else:
	                    regmap = source_dist[min_x:source_dist.shape[0]+1, min_y:source_dist.shape[0]+1]



	                indX2, indY2 = np.where(regmap < searchlensed)


	                reg_centroids = np.zeros( (2, len(indX)))
	                for j in range(len(indX2)):
	                    regmask = np.zeros( (len(regmap[:,0]), len(regmap[0,:])))
	                    region = find_regions(regmap, indX2[j], indY2[j], searchlensed, 0)
	                    regmask[region] = 1.0
	                    reg_centroids[:,j] = centroid(regmap * regmask)


	                #remove duplicates.
	                reg_centroidsx, indexes = np.unique(reg_centroids[0], return_index=True)
	                reg_centroidsy = reg_centroids[1, indexes]


	                n_centroids = 0

	                mindX = np.zeros((len(reg_centroidsx)))
	                mindY = np.zeros((len(reg_centroidsx)))
	                mmu = np.zeros((len(reg_centroidsx)))

	                for j in range(len(reg_centroidsx)):
	                    ic = np.where(np.logical_and(reg_centroids[0] == reg_centroidsx[j], reg_centroids[1] == reg_centroidsy[j]))[0]


	                    if ic.size > 0 and n_centroids == 0:
	                        n_centroids += 1
	#                         print('found multiples')

	                        tempX = np.mean(indX2[ic] + min_x)
	                        tempY = np.mean(indY2[ic] + min_y)
	                        mmu[0] = np.mean(magnification[indX2[ic] + min_x, indY2[ic] + min_y])

	                        mindX[0] = tempX + pos_shift
	                        mindY[0] = tempY + pos_shift

	                    elif ic.size > 0 and n_centroids > 0:
	#                         tempX = mindX
	#                         tempY = mindY
	#                         tempmu = mmu
	#                         n_centroids += 1


	#                         mindX = np.zeros((n_centroids+1))
	#                         mindY = np.zeros((n_centroids+1))
	#                         mmu = np.zeros((n_centroids+1))

	#                         mindX[0:n_centroids] = tempX
	#                         mindY[0:n_centroids] = tempY
	#                         mmu[0:n_centroids] = tempmu

	                        mindX[j] = np.mean(indX2[ic] + min_x) + pos_shift
	                        mindY[j] = np.mean(indY2[ic] + min_y) + pos_shift
	                        m  = np.mean(magnification[indX2[ic] + min_x, indY2[ic] + min_y])


	                        mmu[j] = np.mean(m)

	            elif len(indX) == 1:
	                mindX = np.mean(indX)
	                mindY = np.mean(indY)
	                mmu = np.mean(magnification[int(mindX), int(mindY)])
	                mindX = mindX + pos_shift
	                mindY = mindY + pos_shift
	            else:

	#                 print('None found in SL trying WL')
	                source_x = int(xx[i])
	                source_y = int(yy[i])

	                #these next 4 if statements are to define the position of the WLmask to search for source image
	                #the image will tend to appear outwards the cluster center
	                #the mask is therefore not centered on the source position, but shifted depending on the quadrant
	                if dfxx[i] <= 0 and dfyy[i] <= 0:
	                    xi = int(source_x - WLmask)
	                    xf = int(source_x + WLmask / 2.)
	                    yi 	= int(source_y - WLmask)
	                    yf = int(source_y + WLmask / 2.)
	                elif dfxx[i] < 0 and dfyy[i] > 0:
	                    xi = int(source_x - WLmask / 2.)
	                    xf = int(source_x + WLmask)
	                    yi = int(source_y - WLmask / 2.)
	                    yf = int(source_y + WLmask)
	                elif dfxx[i] > 0 and dfyy[i] > 0:
	                    xi = int(source_x - WLmask / 2.)
	                    xf = int(source_x + WLmask)
	                    yi = int(source_y - WLmask / 2.)
	                    yf = int(source_y + WLmask)
	                elif dfxx[i] > 0 and dfyy[i] < 0:
	                    xi = int(source_x - WLmask / 2.)
	                    xf = int(source_x + WLmask)
	                    yi = int(source_y - WLmask)
	                    yf = int(source_y + WLmask / 2.)
	                #make sure WLmask isnt outside the mask
	                if xi <= 0:
	                    xi = 0
	                if yi <= 0:
	                    yi = 0
	                if xf >= len(wl_dax_dx[:,0])-1:
	                    xf = len(wl_dax_dx[:,0])- 1
	                if yf >= len(wl_dax_dx[0,:])-1:
	                    yf = len(wl_dax_dx[0,:1])-1

	                #not sure if we are concatenating these arrays or adding element by element
	                #estimating magnifications
	                poisson = np.add(wl_dax_dx[xi:xf, yi:yf],wl_day_dy[xi:xf,yi:yf])*scaling_factor
	                wl_daxday_dxdy = np.multiply(wl_dax_dx[xi:xf,yi:yf], wl_day_dy[xi:xf,yi:yf])
	                wl_daxday_dydx = np.multiply(wl_dax_dy[xi:xf,yi:yf], wl_day_dx[xi:xf,yi:yf])
	                magnification = abs(1 / (1 - poisson + np.subtract(wl_daxday_dxdy,wl_daxday_dydx)*scaling_factor ))
	                #find the pixels where the lensed images end up
	                source_dist = np.sqrt( (WLX[xi:xf,yi:yf]-wl_alpha_x[xi:xf,yi:yf]*scaling_factor-source_x)**2.+(WLY[xi:xf,yi:yf]-wl_alpha_y[xi:xf,yi:yf]*scaling_factor-source_y)**2.)
	                indX, indY = np.where(source_dist < searchlensed)

	                if len(indX) > 0:
	                    mmu = np.mean(magnification[indX, indY])
	                    indX += xi
	                    indY += yi

	                    mindX = np.mean(indX)
	                    mindY = np.mean(indY)
	                else:
	                    indX = -999999.0
	                    indY = -999999.0
	                    mu = 0.0
	                    mindX = np.mean(indX)
	                    mindY = np.mean(indY)
	                    mmu = np.mean(mu)

	        else:
	#             print('WL regime')
	            source_x = int(xx[i])
	            source_y = int(yy[i])



	            if dfxx[i] <= 0 and dfyy[i] <= 0:
	                xi = int(source_x - WLmask)
	                xf = int(source_x + WLmask / 2.)
	                yi = int(source_y - WLmask)
	                yf = int(source_y + WLmask / 2.)
	            elif dfxx[i] < 0 and dfyy[i] > 0:
	                xi = int(source_x-WLmask)
	                xf = int(source_x+WLmask/2.)
	                yi = int(source_y-WLmask/2.)
	                yf = int(source_y+WLmask)
	            elif dfxx[i] > 0.0 and dfyy[i] > 0.0:
	                xi = int(source_x-WLmask/2.)
	                xf = int(source_x+WLmask)
	                yi = int(source_y-WLmask/2.)
	                yf = int(source_y+WLmask)
	            elif dfxx[i] > 0.0 and dfyy[i] < 0.0:
	                xi = int(source_x-WLmask/2.)
	                xf = int(source_x+WLmask)
	                yi = int(source_y-WLmask)
	                yf = int(source_y+WLmask/2.)
	            if xi < 0:
	                xi = 0
	            if yi < 0:
	                yi = 0
	            if xf >= len(wl_dax_dx[:,0])-1:
	                xf = len(wl_dax_dx[:,0])-1
	            if yf >= len(wl_dax_dx[0,:])-1:
	                yf = len(wl_dax_dx[0,:])-1

	            xf += 1
	            yf += 1
	            poisson = (wl_dax_dx[xi:xf, yi:yf] + wl_day_dy[xi:xf, yi:yf])*scaling_factor
	            wl_daxday_dxdy = np.multiply(wl_dax_dx[xi:xf,yi:yf], wl_day_dy[xi:xf,yi:yf])
	            wl_daxday_dydx = np.multiply(wl_dax_dy[xi:xf,yi:yf], wl_day_dx[xi:xf,yi:yf])

	            magnification = abs(1 / (1 - poisson + (wl_daxday_dxdy-wl_daxday_dydx)*scaling_factor ))
	            source_dist = np.sqrt( (WLX[xi:xf,yi:yf]-wl_alpha_x[xi:xf,yi:yf]*scaling_factor-source_x)**2.+(WLY[xi:xf,yi:yf]-wl_alpha_y[xi:xf,yi:yf]*scaling_factor-source_y)**2.)

	            indX, indY = np.where(source_dist < searchlensed)
	            if len(indX) > 0:
	                mmu = np.mean(magnification[indX, indY])
	                indX += xi
	                indY += yi

	                mindX = np.mean(indX)
	                mindY = np.mean(indY)
	            else:
	                #if we are here it means there is no image (probably outside of the map)
	                #nonsense position + mag of zero
	                mindX = -999999.0
	                mindY = -999999.0
	                mmu = 0

	        #--- done with WL regime
	        #---- finished lensing analysis


	    if i == 0:
	        x_out = []
	        y_out = []
	        mu_out = []
	        f_out = []


	        if not hasattr(mindX, '__len__'):
	            x_out.append(mindX)
	            y_out.append(mindY)
	            f_out.append(dfflux[i])
	            mu_out.append(mmu)

	        else:
	            for j in range(len(mindX)):
	                x_out.append(mindX[j])
	                y_out.append(mindY[j])
	                mu_out.append(mmu[j])
	                f_out.append(dfflux[i])

	    else:
	        if not hasattr(mindX, '__len__'): #check if an array
	            x_out.append(mindX)
	            y_out.append(mindY)
	            mu_out.append(mmu)
	            f_out.append(dfflux[i])
	        else:
	            for j in range(len(mindX)):
	                x_out.append(mindX[j])
	                y_out.append(mindY[j])
	                mu_out.append(mmu[j])
	                f_out.append(dfflux[i])

	x_out = np.asarray(x_out)
	y_out = np.asarray(y_out)
	f_out = np.asarray(f_out)
	mu_out = np.asarray(mu_out)

	f_final = np.multiply(f_out, mu_out)

	end = time.time()

	print('Finished Lensing sources total time: %s minutes' % str( (end - start) / 60. ) )
	return x_out, y_out, f_final


if __name__ == '__main__':
	catalogue = 'SIDES_PMW_0.npy'
	searchlensed = 0.75
	lens_z = 0.451
	file_dfx = 'alpha_x_ALL_rxj1347_z1p75.txt'
	file_dfy = 'alpha_y_ALL_rxj1347_z1p75.txt'
	LL = 16.6
	lense_x, lense_y, flux = lense(catalogue, lens_z, searchlensed, file_dfx, file_dfy, LL, pixsize=0.25, SLregime=2, main_source_z=1.75, DefField_Norm=None)

	lense_x = lense_x * 0.25 / 6
	lense_y = lense_y * 0.25 / 6

	psf, cf, nc, nbin = get_gaussian_psf_template(3) # assumes pixel fwhm is 3 pixels in each band

	after_lense =  image_model_eval(lense_x, lense_y, nc*flux, 0.0, (166, 166), int(nc), cf)
	plt.title('After Lensing')
	plt.imshow(after_lense, origin='lower')
	plt.colorbar()
	plt.savefig('1600_after_lense.png')
	plt.show()
	plt.clf()

	hdu = fits.PrimaryHDU(after_lense)
	hdul = fits.HDUList([hdu])
	hdul.writeto('full_catalog_after_lense.fits', overwrite=True)
