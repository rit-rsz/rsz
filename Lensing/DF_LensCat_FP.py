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
from read_zitrin import *
import numpy as np
import sys
sys.path.append('../utilities/')
import config
from clus_get_clusparams import *
import matplotlib.pyplot as plt
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


def centroid(arr, start_x, start_y, min, max):
	#not currently being used and not fully tested but saving for later in case.
	ind_x, ind_y = np.where( np.logical_and(arr >= min, arr <= max) == True)

	cent_x = [start_x]
	cent_y = [start_y]

	k = 0
	adj = 1
	bad_sides = 0
	while adj > 0:
		for i in range(cent_x[k] - 1, cent_x[k] + 1):
			for j in range(cent_y[k] - 1, cent_y[k] + 1):
				if i in cent_x and j in cent_y:
					bad_sides += 1

				elif i in ind_x and j in ind_y:
					cent_x.append(i)
					cent_y.append(j)
				else:
					bad_sides += 1
		k += 1
		if k > 12:
			adj = 0 #stop run away program!
		else:
			adj = 8 * len(cent_x) - 8 * k - bad_sides

		print(i)
		print(adj)

	return cent_x, cent_y





def ij(inds, N):
	j = int( inds / n)
	i = round(inds - j * n)
	return i, j


def deriv(y, x):
	#compute the derivative using the definition of a derivative
	# i.e. Delta y / Delta X
	#the derivative function in the IDL code does a 3 point lagrange interpolation and takes the derivative of that
	#it is also weirdly indexed and I'm not sure what the equivalent would be have to play around with it more to
	#triple check that this is giving the same outputs
	dx = np.diff(x) #calculate the difference in x_values
	dy = np.diff(y) #calculate the difference in y_values
	dydx = dy / dx
	return dydx #this function needs to be double checked against the IDL output


def read_input(filename):
	arr  = np.loadtxt(filename, unpack=True)
	return arr

def lenscat(map_size, c_z, pixsize=0.06, searchlensed=0.75, SLregime=[2,2,2], DefField_Norm=None,main_source_z=None,
				file_dfx = 'alpha_x_ALL_rxj1347_z1p75.txt', file_dfy = 'alpha_y_ALL_rxj1347_z1p75.txt',
				mag_file = 'magnification_ALL_rxj1347_z1p75.txt', clusname='rxj1347', Lens_z = 0.4):
	'''
	Inputs:
	map_size : size of the input map
	c_z : the cluster's redshift, Lens_z in Alfredo's code
	pixsize : the size of each pixel in arcsec/pixel, pixscale in Alfredo's code
	searchlensed : default = 0.75 , find the lensed source to extend
	SLregime : Area where we want to apply strong lensing, in arcminN
	'''

	#Parameters to do lensing outside the SL regime --------
	''' We will probably only have one map at a time, so don't need these three...'''
	if SLregime[0] > map_size[0]:
		SLregime[0] == map_size[0]
	if SLregime[1] > map_size[1]:
		SLregime[1] == map_size[1]


	SL_pix = [x * 60. / pixsize for x in SLregime] #converting strong lensing regime from arcmin to pixels
	#for the weakly lensed region we don't estimate the whole image plane but limit to a WLmask x WLmask area around the lensed source
	WLmask = 2.0 # reference scale to estimate WL region in arcmin
	WLmask = WLmask*60./pixsize #WLmask now in pixels

	#IDL version for calculating the luminosity distance / (1 + len_z)**2 @  https://idlastro.gsfc.nasa.gov/ftp/pro/astro/lumdist.pro
	#DL1 = lumdist(Lens_z, H0 = hb, Omega_M = OMm, Lambda0 = OMl,/SILENT) / (1.0+Lens_z)^2, this is the call to lumdist in IDL

	#Python has a version in python @ https://docs.astropy.org/en/stable/cosmology/#functions
	#I assumed a LambdaCDM model based off of what was in Adi and Alfredo's code, but there are other options if that is the wrong model.

	#define a LambdaCDM model with an input hubble constant, matter content, and dark energy content
	cos_model = LCDM(hb*100.0, OMm, OMl)
	# class astropy.cosmology.LambdaCDM(H0, Om0, Ode0, Tcmb0=0, Neff=3.04, m_nu=<Quantity 0. eV>, Ob0=None, name=None)
	'''
	Alternative cos model is directly using built in Planck15
	'''
	DL1 = cos_model.luminosity_distance(c_z) / (1 +c_z)**2

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
		DLS1_DS1 = 1. #if main_source_z was not provided the deflection field normalization should (very likely) be 1


	###########################################################################################################################################


	wl_alpha_x = read_zitrin(map_size,file_dfx)
	wl_alpha_y = read_zitrin(map_size,file_dfy)
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

	YXarray = np.meshgrid(np.arange(wl_lengthy) + 1, np.arange(wl_lengthx) + 1)
	print(YXarray)
	# Not sure which index is needed see below
	WLX = YXarray[0]
	WLY = YXarray[1]

	# Clear some memory
	YXarray = 0

	# this is an example of meshgrid should do
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
	wl_dax_dy = np.array([[np.arange(wl_lengthx)], [np.arange(wl_lengthy)]])
	wl_dax_dx = np.array([[np.arange(wl_lengthx)],[np.arange(wl_lengthy)]])
	wl_day_dy = np.array([[np.arange(wl_lengthx)],[np.arange(wl_lengthy)]])
	wl_day_dx = np.array([[np.arange(wl_lengthx)],[np.arange(wl_lengthy)]])

	# Scipy has something that can maybe do this with misc derivative
	if wl_lengthx == wl_lengthy:
	    for i in range(wl_lengthx):
	        wl_dax_dy = np.diff(wl_alpha_x[i:])
	        wl_dax_dx = np.diff(wl_alpha_x[:i])
	        wl_day_dy = np.diff(wl_alpha_y[i:])
	        wl_day_dx = np.diff (wl_alpha_y[:i])

	else:
	    for i in range(wl_lengthx):
	        wl_dax_dx[:i] = np.diff(wl_alpha_x[:i])
	        wl_day_dx[:i] = np.diff(wl_alpha_y[:i])

	    for i in range(wl_lengthy):
	        wl_dax_dy[i:] = np.diff(wl_alpha_x[i:])
	        wl_day_dy[i:] = np.diff(wl_alpha_y[i:])

	#checking the case in which the arrays are not NxN which we do not have
	# else:
	# 	for i in range(wl_lengthx):
	# 		wl_dax_dx[:i] = diff(wl_alpha_x[:,i], wl_alpha_x[i,:])
	# 		wl_day_dx[:i] = diff(wl_alpha_y[:,i], wl_alpha_y[])
	#
	# 	for i in range(wl_lengthy):
	# 		wl_dax_dy[i:] = diff(wl_alpha_x[i:])
	# 		wl_day_dy[i:] = diff(wl_alpha_y[i:])



# Strong lensing arrays -------------------------------------------------------------------
# here I make a copy of the different arrays (deflection fields, gradients and also the meshgrid arrays)
# but only for the region were we will perform strong lensing analysis. This made thing easier for me when
# I was writting the code and doing tests
	print(wl_lengthx,wl_lengthy)
	'''
	it seems like this whole thing should be done once per band
	Since SL_pix is the size for each band map
	I'm just setting it to 0 for now so I don't have to deal with it
	'''
	print('SL_pix :',SL_pix)

	'''
	wl_lengthx/y is in terms of pixel and SL_pix is in arcseconds, so this won't work
	these two values are always going to make a negative index value
	'''
	dim_minx = int((wl_lengthx-SL_pix[0])/2)
	dim_plusx = int((wl_lengthx+SL_pix[0])/2)
	dim_miny = int((wl_lengthy-SL_pix[0])/2)
	dim_plusy = int((wl_lengthy+SL_pix[0])/2)
	print(dim_minx, dim_plusx, dim_miny, dim_plusy)

	print(wl_lengthx, wl_lengthy)
	print(SL_pix)
	print(wl_alpha_x.shape)
	# This seems like it might not be the right python syntax
	sl_alpha_x = wl_alpha_x[dim_minx:dim_plusx,dim_minx:dim_plusx]
	sl_alpha_y = wl_alpha_y[dim_miny:dim_plusy,dim_miny:dim_plusy]

	print(wl_dax_dy.shape)
	sl_dax_dy = wl_dax_dy[dim_miny:dim_plusy,dim_miny:dim_plusy]
	sl_dax_dx = wl_dax_dx[dim_minx:dim_plusx,dim_minx:dim_plusx]
	sl_day_dy = wl_day_dy[dim_miny:dim_plusy,dim_miny:dim_plusy]
	sl_day_dx = wl_day_dx[dim_minx:dim_plusx,dim_minx:dim_plusx]
	sl_lengthx = len(sl_alpha_x[:0])

	YXarray = np.meshgrid(np.arange(sl_lengthx)+1, np.arange(sl_lengthx)+1)
	SLX = YXarray[0]
	SLY = YXarray[1]
	YXarray = 0

	# shift between WL and SL frames (you will need this number to estimate the correct positions of the weak and strong lensed sources)
	pos_shift = (wl_lengthx-sl_lengthx)/2.

	#------------------------------------------------------------- OJO -----------------------------------------------------------------

	# Read source catalogue (input_cat)
	read_default_format = 'l,d,d,f,f'   #Default Format
	print('Reading Catalogue to be lensed')
	# Uses read call but this kind of strange. This is I think reading lines as varriables like dict keys
	# readcol, input_cat, dfid, dfxx, dfyy, dfredz, dfflux, skipline = 1, FORMAT = read_default_format
	filename = config.CLUSDATA + 'sides_sims/' + clusname + '/SIDES_PLW_sim1.npy'
	truth_table = np.load(filename ,allow_pickle=True)

	dfxx = truth_table.item().get('RA')
	dfyy = truth_table.item().get('DEC')
	dfflux = truth_table.item().get('Flux')
	dfredz = truth_table.item().get('Redshift')
	# Convert coordinates from arcsecs to pixels in the larger WL frame
	xx = (dfxx + map_size[0] * 60./2.) / pixsize
	yy = (dfyy + map_size[1] * 60./2.) / pixsize
	# Number of sources in catalogue
	nsources = len(xx)


	for i in range(nsources):
		# begin dealing with foreground sources (we don't need to lens these) ---------------------
		if Lens_z > dfredz[i]:
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
			DLS2 = (DS2*(1.+dfredz[i]) - DL1*(1.+Lens_z))
			#this is the factor to scale the deflection field to different source redshifts.
			scaling_factor = (DLS2/DS2) / DLS1_DS1
			# ; begin SL regime analysis

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
						if (min_x + multi_ll) < len(source_dist[:,0]):
							if min_y + multi_ll < len(source_dist[0,:]):
								regmap = source_dist[min_x:min_x + multi_ll, min_y:min_y+multi_ll]
							else:
								regmap = source_dist[min_x:min_x + multi_ll, min_y-1:len(source_dist[:,0])-1]
						elif min_y + multi < len(source_dist[0,:]):
							regmap = source_dist[min_x-1:len(source_dist[:,0]), min_y:min_y + multi_ll]
						else:
							regmap = source_dist[min_x-1:len(source_dist[:,0])-1, min_y-1:len(source_dist[0,:])-1]
						indXY2 = np.where(regmap < searchlensed)
						indXY2 = ij( indXY2, len(regmap[:,0]))
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
							#find the centroids within the masked region
							reg_centroids[:, j] = placeholder#need a function to do this can't find equiv
						n_centroids = 0 #this is a counter initialized at 0.
						for j in range(len(reg_centroids[0,:])):
							i1 = np.where(reg_centroids[0,:] == reg_centroids[0,j])
							i2 = np.where(reg_centroids[1,:] == reg_centroids[1,j])
							ic = np.intersect1d(i1, i2) #find whree i1 and i2 overlap (effectivley an and operator)
							if len(ic) != 0 and n_centroids > 0:
								print('we found multiples')
								n_centroids += 1
								j = j + len(ic) - 1 #dont double count centroids
								#sav mindX, mindY and magnification of each image
								temp_indX = mindX
								temp_indY = mindY
								temp_mu = mmu
								mindX = np.zeros(len(n_centroids))
								mindY = np.zeros(len(n_centroids))
								mmu = np.zeros(len(n_centroids))
								mindX[0:n_centroids-2] = temp_indX
								mindX[n_centroids-1] = np.mean(indXY2[ic,0] + min(indXY[:,0])) + pos_shift
								mindY[0:n_centroids-2] = temp_indY
								mindY[n_centroids-1] = np.mean(indXY2[ic,1] + min(indXY[:,1])) + pos_shift
								mmu[0:n_centroids-2] = temp_mu
								mmu[n_centroids-1] = np.mean(magnification[indXY2[ic, 0] + min(indXY[:,0]), indXY2[ic, 1] + min(indXY[:,1])])
							elif ic[0] != -1 and n_centroids == 0:
								n_centroids += 1
								j = j + len(ic) - 1 #no double counting centroids !
								mindX = np.mean(indXY2[ic,0] + min(indXY[:,0]))
								mindY = np.mean(indXY2[ic,1] + min(indXY[:,1]))
								#on line 351 (322 in IDL) magnification was being indexed but below here it was being called as if it was
								#a function.
								mmu = np.mean(magnification[indXY2[ic, 0] + min(indXY[:,0]), indXY2[ic, 1] + min(indXY[:,1])])
								mindX = mindX + pos_shift
								mindY = mindY + pos_shift
					else:
						mindX = np.mean(indXY[:,0])
						mindY = np.mean(indXY[:,1])
						mmu = np.mean(magnification[mindX, mindY])
						mindX = mindX + pos_shift
						mindY = mindY + pos_shift
				else:
					print('No images found in the SL regime, now checking for WL image')
					source_x = xx[i]
					source_y = yy[i]
					if dfxx[i] <= 0 and dfyy[i] <= 0:
						xi = source_x - WLmask
						xf = source_x + WLmask / 2.
						yi 	= source_y - WLmask
						yf = source_y + WLmask / 2.
					elif dfxx[i] < 0 and dfyy[i] > 0:
						xi = source_x - WLmask / 2.
						xf = source_x + WLmask
						yi = source_y - WLmask / 2.
						yf = source_y + WLmask
					elif dfxx[i] > 0 and dfyy[i] > 0:
						xi = source_x - WLmask / 2.
						xf = source_x + WLmask
						yi = source_y - WLmask / 2.
						yf = source_y + WLmask
					elif dfxx[i] > 0 and dfyy[i] < 0:
						xi = source_x - WLmask / 2.
						xf = source_X + WLmask
						yi = source_y - WLmask
						yf = source_y + Wlmask / 2.
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
					poisson = (wl_dax_dx[xi:xf, yi:yf] + wl_day_dy[xi:xf,yi:yf])*scaling_factor
					magnification = abs(1 / (1 - poisson + (wl_dax_dx[xi:xf,yi:yf]*wl_day_dy[xi:xf,yi:yf]-wl_dax_dy[xi:xf,yi:yf]*wl_day_dx[xi:xf,yi:yf])*scaling_factor ))
					#find the pixels where the lensed images end up
					source_dist = sqrt( (WLX[xi:xf, yi:yf]-wl_alpha_x[xi:xf, yi:yf]-wl_alpha_x[xi:xf, yi:yf]*scaling_factor-source_x)**2 +(WLY[xi:xf,yi:yf]-wl_alpha_y[xi:xf,yi:yf]*scaling_factor-source_y)**2)
					indXY = np.where(source_dist < searchlensed)
					if indXY[0] != -1:
						indXY = ij(indXy, len(poisson[:,0])) #still need to find a function for ij
						indX = indXY[:,0]
						indY = indXY[:,1]
						mu = magnification(indX, indY) #so magnification is also a function?
						indX = indXY[:,0] + xi
						indY = indXY[:,1] + yi
					else:
						indX = -999999.0
						indY = -999999.0
						mu = 0.0
					mindX = np.mean(indX)
					mindY = np.mean(indY)
					mmu = np.mean(mu)
			else:
				print('WL regime')
				source_x = xx[i]
				source_y = yy[i]
				if dfxx[i] <= 0 and dfyy[i] <= 0:
					xi = source_x - WLmask
					xf = source_X + WLmask / 2.
					yi = source_y - WLmask
					yf = source_y + WLmask / 2.
				elif dfxx[i] < 0 and dfyy[i] > 0:
					xi = source_x-WLmask
					xf = source_x+WLmask/2.
					yi = source_y-WLmask/2.
					yf = source_y+WLmask
				elif dfxx[i] > 0.0 and dfyy[i] > 0.0:
					xi = source_x-WLmask/2.
					xf = source_x+WLmask
					yi = source_y-WLmask/2.
					yf = source_y+WLmask
				elif dfxx[i] > 0.0 and dfyy[i] < 0.0:
					xi = source_x-WLmask/2.
					xf = source_x+WLmask
					yi = source_y-WLmask
					yf = source_y+WLmask/2.
				if xi < 0:
					xi = 0
				if yi < 0:
					yi = 0
				if xf >= len(wl_dax_dx[:,0])-1:
					xf = len(wl_dax_dx[:,0])-1
				if yf >= len(wl_dax_dx[0,:])-1:
					yf = len(wl_dax_dx[0,:])-1
				poisson = (wl_dax_dx[xi:xf, yi:yf] + wl_day_dy[xi:xf, yi:yf])*scaling_factor
				magnification = abs(1 / (1 - poisson + (wl_dax_dx[xi:xf,yi:yf]*wl_day_dy[xi:xf,yi:yf]-wl_dax_dy[xi:xf,yi:yf]*wl_day_dx[xi:xf,yi:yf])*scaling_factor ))
				source_dist = sqrt( (WLX[xi:xf,yi:yf]-wl_alpha_x[xi:xf,yi:yf]*scaling_factor-source_x)**2.+(WLY[xi:xf,yi:yf]-wl_alpha_y[xi:xf,yi:yf]*scaling_factor-source_y)**2.)
				indXY = np.where(source_dist < searchlensed)
				if indXY[0] != -1:
					indXY = ij(indXy, len(poisson)) #still need this function equivalent
					indX = indXY[:,0]
					indY = indXY[:,1]
					mu = magnification(indX, indY)
					indX = indXy[:,0] + xi
					indY = indXY[:,1] + yi
				else:
					#if we are here it means there is no image (probably outside of the map)
					#nonsense position + mag of zero
					indX = -999999.0
					indY = -999999.0
					mu = 0
					mindX = np.mean(indX)
					mindY = np.mean(indY)
					mmu = np.mean(mu)
					#--- done with WL regime
					#---- finished lensing analysis
		if i == 0:
	#constructing output arrays
			mindX = int(mindX)
			mindY = int(mindY)
			x_out = np.arange(mindX)
			y_out = np.arange(mindY)
			mu_out = mmu
			x_in = np.zeros(mindX)
			x_in[:] = dfxx[i]
			y_in = np.zeros(mindX)
			y_in[:] = dfyy[i]
			z_in = np.zeros(mindX)
			z_in[:] = dfredz[i]
			f_in = np.zeros(mindX)
			f_in[:] = dfflux[i]
		else:

			mindX = int(mindX)
			mindY = int(mindY)

			x_temp = np.zeros( mindX + len(x_out))
			x_temp[0:len(x_out)-1] = len(x_out)
			x_temp[len(x_out):-1] = mindX
			x_out = x_temp
			y_temp = np.zeros(mindY + len(y_out))
			y_temp[0:len(y_out)-1] = len(y_out)
			y_temp[len(y_out):-1] = mindY
			y_out = y_temp

			xi_temp = np.zeros(mindX + len(x_in))
			xi_temp[0:len(x_in)] = x_in
			xi_temp[len(x_in):-1] = dfxx[i]
			x_in = xi_temp
			yi_temp = np.zeros(mindX + len(y_in))
			yi_temp[0:len(y_in)] = y_in
			yi_temp[len(y_in):-1] = dfyy[i]
			y_in = yi_temp
			fi_temp = np.zeros(mindX + len(f_in))
			fi_temp[0:len(f_in)] = f_in
			fi_temp[len(f_in):-1] = dfflux[i]
			f_in = fi_temp
	print('done')


	return x_in,y_in, z_in,f_in, x_out ,y_out, mu_out
if __name__ == '__main__':

	params, err = clus_get_clusparams('rxj1347')
	# cat_file =  '/data/vaughan/SPIRE/lensing_test_catalogues/npy_cats/SIDES_PSW_sim0.npy'
	cat = np.load(config.HOME + 'sides_sims/sides_PSW_sim1.npy' ,allow_pickle=True)
	# cat = np.load(cat_file, allow_pickle=True)
	#map_size = cat.item().get('RA').shape this won't work need to find mapsize diff way
	map_size = [290,270]
	c_z = params['z']
	lenscat(map_size, c_z, pixsize=6, searchlensed=0.75, SLregime=[2,2], DefField_Norm=None,main_source_z=None,
					file_dfx = 'IDL_program/alpha_x_ALL_rxj1347_z1p75.txt', file_dfy = 'IDL_program/alpha_y_ALL_rxj1347_z1p75.txt',
					mag_file = 'magnification_ALL_rxj1347_z1p75.txt')
