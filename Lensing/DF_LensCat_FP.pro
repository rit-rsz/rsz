;+
; NAME: DF_LensCat_FP
;
; PURPOSE: function called by DF_MakeMaps in order to lens a catalogue using deflection fields
;
;
; INPUTS:	input_cat: input catalogue to be lensed.
;			FILE_DFX / FILE_DFY: deflection field in X & Y (provided by Adi Zitrin. Use READ_ZITRIN.pro read this files).
;			Lens_z: cluster (lens) redshift.
;			Pix_Scale: pixel size (Pix_Scale x Pix_Scale) in arcsec (default: 0.06" x 0.06").
;			LL: size of square map (L x L) in arcmin.
;			NPIXELS: number of pixels in map (NPIXELS).
;			SEARCHLENSED: this number is to find the lensed source which can be extended (arcs). Use the default value 0.75
;
; OPTIONAL INPUTS:  SLREGIME [arcmin] for scales larger than SLREGIME we will not estimate the entire image plane. This makes the code faster.
;							Default value is 2 arcmin.
;					MAIN_SOURCE_Z (optional): redshift of source used to estimate deflections fields (provided by Adi Zitrin). If not provided
;							the code assumes the lensing normalization is 1 or equal to DEFFIELD_NORM (if provided)
;					DEFFIELD_NORM (optional): deflection field normalization (provided by Adi Zitrin?). If not provided (and MAIN_SOURCE_Z was)
;							also not provided), the normalization is assumed to be 1.
;
; OUTPUTS: an array with [source ID (input), X position (input), Y position (input), redshift (input),
;							unlensed flux (input), lensed X position, lensed Y position, magnification]
;			Note that, to build the lensed map, you have to you should use the lensed positions and unlensed flux*magnification
;
;  Created by AMB on 2018-05-08.
;-

@meshgrid.pro
@READ_ZITRIN.pro

FUNCTION DF_LensCat_FP, input_cat, FILE_DFX, FILE_DFY, Lens_z, Pix_Scale, LL, NPIXELS=npixels,MAIN_SOURCE_Z = main_source_z, DEFFIELD_NORM = DefField_Norm, SEARCHLENSED = searchlensed, SLREGIME = SLRegime


;COSMOLOGY

	hb  = 0.7                    ; Hubble constante (Ho = 100*hb)
	OMm = 0.3                    ; Omega Matter (dark + baryonic)
	OMl = 0.7                    ; Omega Lambda
	OMk = 1. - OMm - OMl          ; Omega curvature
	DH = 2997.92/hb                ; Hubble Distance in Mpc => DH = (c * 1e-3)/Ho

;	;using Planck 2015 cosmology ... Ade+2015 a (Table 4 last column)
;	hb  = 0.6774                    ; Hubble constante (Ho = 100*hb)
;	OMm = 0.3089                    ; Omega Matter (dark + baryonic)
;	OMl = 0.6911                    ; Omega Lambda
;	OMk = 1.0 - OMm - OMl           ; Omega curvature
;	DH  = 2997.92/hb                  ; Hubble Distance in Mpc => DH = (c * 1e-3)/Ho

;Define input cluster redshift, source redshift to which deflection field is scaled to, and the pixel
;scale (both are input provided by Adi Zitrin). Often the deflection field will be scaled to DLS/DS=1.

;if pix_scale is not provided we use the default value of 0.06 arcsec
IF n_elements(pix_scale) LT 1 THEN BEGIN
	pix_scale=0.06          ;arcsec/pixel
ENDIF

;if searchlensed is not provided we use the default value of 0.75
IF n_elements(searchlensed) LT 1 THEN BEGIN
	searchlensed=0.75
ENDIF

;Parameters to do lensing outside the SL regime --------

;if SLregime is not provided we will assume that Strong Lensing is within the central 2 x 2 s    #DL1 = lumdist(Lens_z, H0 = hb, Omega_M = OMm, Lambda0 = OMl,/SILENT) / (1.0+Lens_z)^2
q.arcmin
IF n_elements(SLregime) LT 1 THEN SLregime = 2.0 ;arcmin

;if SLregime is larger than the map size (LL) we will set SLregime = LL
IF SLregime GT LL THEN SLregime = LL
SLregime_pix = SLregime*60./pix_scale  ;SL regime now in pixels

;for the weakly lensed region we don't estimate the whole image plane but limit to a WLmask x WLmask area around the lensed source
WLmask = 2.0 ; reference scale to estimate WL region in arcmin
WLmask = WLmask*60./pix_scale ; WLmask now in pixels


;estimate angular diameter distance to the lensing cluster (z = Lens_z)
;lumdist is an IDL function that returns the luminosity distance (DL) in Mpc. Angular diameter distance DA = DL/(1+z)^2
DL1 = lumdist(Lens_z, H0 = hb, Omega_M = OMm, Lambda0 = OMl,/SILENT) / (1.0+Lens_z)^2


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

;if instead of providing main_source_z, a normalization (DefField_Norm) different to 1 was provided here we assign it to the varibale DLS1_DS1
IF n_elements(DefField_Norm) EQ 1 THEN DLS1_DS1 = DefField_Norm

;we read the lensing deflection fields provided by Adi Zitrin, and which have npixels x npixels elements as a text file
;we use the function READ_ZITRIN which takes as input the number of pixels in one side of the map (npixels) and the name of the file
	wl_alpha_x = READ_ZITRIN(npixels,FILE_DFX)     ;read the deflection field in X
	wl_alpha_y = READ_ZITRIN(npixels,FILE_DFY)     ;read the deflection field in Y

	wl_alpha_x = TRANSPOSE(wl_alpha_x)             ;need to transpose (at least in IDL)
	wl_alpha_y = TRANSPOSE(wl_alpha_y)             ;need to transpose (at least in IDL)

;make sure that the deflections fields are NxN
IF N_ELEMENTS(wl_alpha_x[*,0]) GE N_ELEMENTS(wl_alpha_x[0,*]) THEN BEGIN
	wl_alpha_x = wl_alpha_x[0:N_ELEMENTS(wl_alpha_x[0,*])-1,*]
	wl_alpha_y = wl_alpha_y[0:N_ELEMENTS(wl_alpha_y[0,*])-1,*]
ENDIF

IF N_ELEMENTS(wl_alpha_x[*,0]) LT N_ELEMENTS(wl_alpha_x[0,*]) THEN BEGIN
	wl_alpha_x = wl_alpha_x[*,0:N_ELEMENTS(wl_alpha_x[*,0])-1]
	wl_alpha_y = wl_alpha_y[*,0:N_ELEMENTS(wl_alpha_y[*,0])-1]
ENDIF

;Weak lensing: estimate X & Y arrays
wl_lengthx = N_ELEMENTS(wl_alpha_x[*,0])
wl_lengthy = N_ELEMENTS(wl_alpha_x[0,*])
PRINT, 'This should be 1 ... ', wl_lengthx/wl_lengthy	;just to check

YXarray=meshgrid(FINDGEN(wl_lengthy)+1.0,FINDGEN(wl_lengthx)+1.0)  ;this function is equivalent to meshgrid(x,y) in matlab and I believe NumPy also has it too
WLX = YXarray[*,*,0]
WLY = YXarray[*,*,1]
YXarray = 0     ; get rid of this array just to save memory

;this is an example of what meshgrid does
;a =  MESHGRID(findgen(3),findgen(3))
;print, a[*,*,0]
;      0.00000      1.00000      2.00000
;      0.00000      1.00000      2.00000
;      0.00000      1.00000      2.00000
;print, a[*,*,1]
;      0.00000      0.00000      0.00000
;      1.00000      1.00000      1.00000
;      2.00000      2.00000      2.00000

;Weak lensing: estimate gradients
wl_dax_dy = FLTARR(wl_lengthx,wl_lengthy)      ;generate empty arrays were we'll save the derivatives
wl_dax_dx = FLTARR(wl_lengthx,wl_lengthy)      ;generate empty arrays were we'll save the derivatives
wl_day_dy = FLTARR(wl_lengthx,wl_lengthy)      ;generate empty arrays were we'll save the derivatives
wl_day_dx = FLTARR(wl_lengthx,wl_lengthy)      ;generate empty arrays were we'll save the derivatives

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

;Strong lensing arrays -------------------------------------------------------------------
;here I make a copy of the different arrays (deflection fields, gradients and also the meshgrid arrays)
;but only for the region were we will perform strong lensing analysis. This made thing easier for me when
;I was writting the code and doing tests

	sl_alpha_x = wl_alpha_x[(wl_lengthx-SLregime_pix)/2.:(wl_lengthx+SLregime_pix)/2.,(wl_lengthx-SLregime_pix)/2.:(wl_lengthx+SLregime_pix)/2.]
	sl_alpha_y = wl_alpha_y[(wl_lengthx-SLregime_pix)/2.:(wl_lengthx+SLregime_pix)/2.,(wl_lengthx-SLregime_pix)/2.:(wl_lengthx+SLregime_pix)/2.]

	sl_dax_dy = wl_dax_dy[(wl_lengthx-SLregime_pix)/2.:(wl_lengthx+SLregime_pix)/2.,(wl_lengthx-SLregime_pix)/2.:(wl_lengthx+SLregime_pix)/2.]
	sl_dax_dx = wl_dax_dx[(wl_lengthx-SLregime_pix)/2.:(wl_lengthx+SLregime_pix)/2.,(wl_lengthx-SLregime_pix)/2.:(wl_lengthx+SLregime_pix)/2.]
	sl_day_dy = wl_day_dy[(wl_lengthx-SLregime_pix)/2.:(wl_lengthx+SLregime_pix)/2.,(wl_lengthx-SLregime_pix)/2.:(wl_lengthx+SLregime_pix)/2.]
	sl_day_dx = wl_day_dx[(wl_lengthx-SLregime_pix)/2.:(wl_lengthx+SLregime_pix)/2.,(wl_lengthx-SLregime_pix)/2.:(wl_lengthx+SLregime_pix)/2.]

	sl_lengthx = N_ELEMENTS(sl_alpha_x[*,0])

	YXarray=meshgrid(FINDGEN(sl_lengthx)+1.0,FINDGEN(sl_lengthx)+1.0)
	SLX = YXarray[*,*,0]
	SLY = YXarray[*,*,1]
	YXarray = 0

;shift between WL and SL frames (you will need this number to estimate the correct positions of the weak and strong lensed sources)
pos_shift = (wl_lengthx-sl_lengthx)/2.

;------------------------------------------------------------- OJO -----------------------------------------------------------------

;Read source catalogue (input_cat)
read_default_format = 'l,d,d,f,f'      ;Default format
print, 'Reading Catalogue to be lensed'
readcol, input_cat, dfid, dfxx, dfyy, dfredz, dfflux, skipline = 1, FORMAT = read_default_format
;dfid=source ID / dfxx=x position in arcsec / dfyy=y position in arcsec / dfredz=source redshift / dfflux=source flux

	;convert coordinates from arcsecs to pixels in the larger WL frame
	xx = LONG((dfxx+LL*60./2.)/Pix_Scale) ;+ 1
	yy = LONG((dfyy+LL*60./2.)/Pix_Scale) ;+ 1
	;Number of sources in catalogue
	NSources = N_ELEMENTS(xx)

	FOR ii=0.,NSources-1. DO BEGIN
		ti = systime(/seconds)          ; don't worry about this, it's only to know when the lensing analysis started
		PRINT, 'Lensing source = ', ii+1, ' / ', NSources      ; don't worry about this, it just prints which source you are dealing with

; begin dealing with foreground sources (we don't need to lens these) ---------------------
		IF (Lens_z GE dfredz[ii]) THEN BEGIN
			print, 'Source ', ii+1,' is a foreground source.'
			mindX = xx[ii]
			mindY = yy[ii]
			mmu = 1.0
		ENDIF ELSE BEGIN
; finish dealing with foreground sources --------------------------------------------------

; begin dealing with background sources (we do need to lens these)
			;estimate angular diameter distance to source (z = dfredz)
			DS2 = lumdist(dfredz[ii], H0 = hb, Omega_M = OMm, Lambda0 = OMl,/SILENT) / (1.0+dfredz[ii])^2.
			;estimate distance between the lensing cluster and source
			DLS2 = (DS2*(1.+dfredz[ii]) - DL1*(1.+Lens_z))/(1.+dfredz[ii])     ;WARNING: this is only valid for Omega_k = 0

			;this is the factor to scale the deflection field to different source redshifts.
			scaling_factor = (DLS2/DS2) / DLS1_DS1

; begin SL regime analysis
			print, 'time 1 = ', (systime(/seconds) - ti)/60.               ;don't worry about this
;			IF SQRT(dfxx[ii]^2.+dfyy[ii]^2.) LE SLregime*60. THEN BEGIN    ;don't worry about this

			;check if source is within strong lensing region
			IF ( ABS(dfxx[ii]) LE SLregime*60./2. ) AND ( ABS(dfyy[ii]) LE SLregime*60./2. ) THEN BEGIN
			print, 'SL regime'

			;estimate source position in the SL frame
			source_x = LONG((dfxx[ii]+SLregime*60./2.)/Pix_Scale) ;might need to add + 1
			source_y = LONG((dfyy[ii]+SLregime*60./2.)/Pix_Scale) ;might need to add + 1

			;estimate magnifications (these come from lensing equations)
			poisson = DOUBLE(sl_dax_dx+sl_day_dy)*scaling_factor
			magnification=ABS(1./(1.-poisson+DOUBLE(sl_dax_dx*sl_day_dy-sl_dax_dy*sl_day_dx)*scaling_factor))

			;find pixels where the lensed image (or images) end up
			source_dist=sqrt(DOUBLE(SLX-sl_alpha_x*scaling_factor-source_x)^2.+DOUBLE(SLY-sl_alpha_y*scaling_factor-source_y)^2.)
			indXY = WHERE(source_dist LT searchlensed)


			;if a lensed source was found (within searchlensed) ...
			IF indXY[0] NE -1 THEN BEGIN
				indXY = ij(DOUBLE(indXY), sl_lengthx)          ; find the (i,j) indices corresponding to the image pixels.
print,' -----------------1 ', indXY				;don't worry about this
;here we check for multiplicity ---------------------------------------------------------------------
;if more than 1 pixel satisfy "source_dist LT searchlensed" we might have multiple images of a source
				IF N_ELEMENTS(indXY[*,0]) GT 1 THEN BEGIN

				;we cut a square sub-map including all pixels with "source_dist LT searchlensed", so we don't have to deal with the whole map
				temp_multi_x = max(indXY[*,0]) - min(indXY[*,0])
				temp_multi_y = max(indXY[*,1]) - min(indXY[*,1])

				IF (temp_multi_x - temp_multi_y) GE 0.0 THEN multi_ll = temp_multi_x ELSE multi_ll = temp_multi_y

				;a bunch of IFs to consider cases where the lensed source lies close to the edge of the map
				IF (min(indXY[*,0])+multi_ll) LT N_ELEMENTS(source_dist[*,0]) THEN BEGIN
					IF (min(indXY[*,1])+multi_ll) LT N_ELEMENTS(source_dist[0,*]) THEN BEGIN
						regmap = source_dist[min(indXY[*,0]):min(indXY[*,0])+multi_ll,min(indXY[*,1]):min(indXY[*,1])+multi_ll]
					ENDIF ELSE BEGIN
						regmap = source_dist[min(indXY[*,0]):min(indXY[*,0])+multi_ll,min(indXY[*,1])-1.:N_ELEMENTS(source_dist[*,0])-1.]
					ENDELSE
				ENDIF ELSE BEGIN
					IF (min(indXY[*,1])+multi_ll) LT N_ELEMENTS(source_dist[0,*]) THEN BEGIN
						regmap = source_dist[min(indXY[*,0])-1.:N_ELEMENTS(source_dist[*,0])-1.,min(indXY[*,1]):min(indXY[*,1])+multi_ll]
					ENDIF ELSE BEGIN
						regmap = source_dist[min(indXY[*,0])-1.:N_ELEMENTS(source_dist[*,0])-1.,min(indXY[*,1])-1.:N_ELEMENTS(source_dist[*,0])-1.]
					ENDELSE
				ENDELSE

				indXY2 = where(regmap LT searchlensed)                  ;find pixels within the sub-map where the lensed image (or images) end up
print,' -----------------2 ', indXY2									;don't worry about this
				indXY2 = ij(DOUBLE(indXY2), N_ELEMENTS(regmap[*,0]))	;find the (i,j) indices corresponding to the image pixels.

				;here we look for contiguous regions (within sub-map) that satisfy "source_dist LT searchlensed"
				reg_centroids = FLTARR(2,N_ELEMENTS(indXY[*,0]))      ;in this array we will save the (i,j) pairs

				FOR jj=0., N_ELEMENTS(indXY2[*,0])-1. DO BEGIN
					;find all pixels contiguous to pixel indXY2[jj,0], indXY2[jj,1] that satisfy searchlensed condition (IDL routine)
					region = SEARCH2D(regmap, indXY2[jj,0], indXY2[jj,1], 0, searchlensed, /DIAGONAL)
					;make an array the size of the sub-map with all pixels=0 to construct a mask
					regmask = FLTARR(N_ELEMENTS(regmap[*,0]),N_ELEMENTS(regmap[*,0]))
					;give a value of 1 to all pixels contiguous to pixel indXY2[jj,0], indXY2[jj,1] that satisfy searchlensed condition
					regmask[region] = 1.0
					;find the centroid of the contiguous region and save its (i,j) values in reg_centroids
					reg_centroids[*,jj] = centroid(regmap*regmask)
					;note that reg_centroids will have repeated centroids. We clean this array below.
				ENDFOR

				n_centroids = 0.0            ;define a counter where we'll save how many images a source has

				FOR jj=0., N_ELEMENTS(reg_centroids[0,*])-1. DO BEGIN
					;we check within reg_centroids how many times we've repeated the same centroid
					ic = where(reg_centroids[0,*] EQ reg_centroids[0,jj] AND reg_centroids[1,*] EQ reg_centroids[1,jj], COMPLEMENT = ic2)

					IF (ic[0] NE -1) AND (n_centroids GT 0.0) THEN BEGIN
					;IF ic = -1 means that we only have 1 image of this source (see IF below)
					;IF n_centroid is GT 0 it means we've already find 1 image (see IF below) and we have multiple images
						print, 'we found multiples'
						n_centroids = n_centroids+1.0       ;to keep track of how many times we've entered this loop
						jj = jj + N_ELEMENTS(ic) - 1.0      ;skip all other elements of reg_centroids that correspond to this same centroid


						;in the following lines we save the position (mindX,mindY) and magnification (mmu) of each image.
						temp_indX = mindX
						temp_indY = mindY
						temp_mu   = mmu

						mindX = fltarr(n_centroids)
						mindY = fltarr(n_centroids)
						mmu   = fltarr(n_centroids)

						mindX[0:n_centroids-2.] = temp_indX
						mindX[n_centroids-1.] = MEAN(indXY2[ic,0]+min(indXY[*,0])) + pos_shift

						mindY[0:n_centroids-2.] = temp_indY
						mindY[n_centroids-1.] = MEAN(indXY2[ic,1]+min(indXY[*,1])) + pos_shift

						mmu[0:n_centroids-2.] = temp_mu
						mmu[n_centroids-1.] = MEAN( magnification[ indXY2[ic,0]+min(indXY[*,0]),indXY2[ic,1]+min(indXY[*,1]) ] )
					ENDIF

					IF (ic[0] NE -1) AND (n_centroids EQ 0.0) THEN BEGIN
					;IF ic = -1 there is something wrong going on (so this condition is just to check)
					;IF n_centroid is EQ 0 it means we have at least 1 image of this source, so we add 1 to the centroid counter (n_centroids)
						n_centroids = n_centroids+1.0
						jj = jj + N_ELEMENTS(ic) - 1.			;skip all other elements of reg_centroids that correspond to this same centroid

						;in the following lines we save the position (mindX,mindY) and magnification (mmu) of each image.
						mindX = MEAN(indXY2[ic,0]+min(indXY[*,0]))
						mindY = MEAN(indXY2[ic,1]+min(indXY[*,1]))
						mmu   = MEAN( magnification( indXY2[ic,0]+min(indXY[*,0]),indXY2[ic,1]+min(indXY[*,1]) ) )
						mindX = mindX + pos_shift
						mindY = mindY + pos_shift
					ENDIF

				ENDFOR

				ENDIF ELSE BEGIN ; this should be the case when only one pixel of source_dist LT searchlensed (i.e. only one image)
					mindX = MEAN(indXY[*,0])
					mindY = MEAN(indXY[*,1])
					mmu   = MEAN(magnification[mindX,mindY])
					mindX = mindX + pos_shift
					mindY = mindY + pos_shift
				ENDELSE

;finished checking for multiplicity ---------------------------------------------------------------------
;If we do not find an image for a source in the SL regime (e.g. it is close to the egde between SL and WL maps), then we do WL analysis
				ENDIF ELSE BEGIN
					print, 'No image found in the SL regime ... checking for WL image.'         ;don't worry about this
					source_x = xx[ii]       ;this is just the X position of the source in pixel units
					source_y = yy[ii]       ;this is just the Y position of the source in pixel units

					;The next 4 IFs are to define the position of the WLmask to search for the source image.
					;The image will tend to appear outwards the cluter center (at least in the WL regime).
					;The mask is therefore not centered on the source position, but shifted depending on the quadrant where the source is.
					IF dfxx[ii] LE 0.0 AND dfyy[ii] LE 0.0 THEN BEGIN
						xi = source_x-WLmask
						xf = source_x+WLmask/2.
						yi = source_y-WLmask
						yf = source_y+WLmask/2.
					ENDIF

					IF dfxx[ii] LT 0.0 AND dfyy[ii] GT 0.0 THEN BEGIN
						xi = source_x-WLmask
						xf = source_x+WLmask/2.
						yi = source_y-WLmask/2.
						yf = source_y+WLmask
					ENDIF

					IF dfxx[ii] GT 0.0 AND dfyy[ii] GT 0.0 THEN BEGIN
						xi = source_x-WLmask/2.
						xf = source_x+WLmask
						yi = source_y-WLmask/2.
						yf = source_y+WLmask
					ENDIF

					IF dfxx[ii] GT 0.0 AND dfyy[ii] LT 0.0 THEN BEGIN
						xi = source_x-WLmask/2.
						xf = source_x+WLmask
						yi = source_y-WLmask
						yf = source_y+WLmask/2.
					ENDIF

					;These 4 IFs are just to make sure the WLmask do not fall outside the map
					IF xi LE 0.0 THEN xi = 0.0
					IF yi LE 0.0 THEN yi = 0.0
					IF xf GE N_ELEMENTS(wl_dax_dx[*,0])-1. THEN xf = N_ELEMENTS(wl_dax_dx[*,0])-1.
					IF yf GE N_ELEMENTS(wl_dax_dx[0,*])-1. THEN yf = N_ELEMENTS(wl_dax_dx[0,*])-1.

				;estimate magnifications (these come from lensing equations)
				poisson = DOUBLE(wl_dax_dx[xi:xf,yi:yf] + wl_day_dy[xi:xf,yi:yf])*scaling_factor
				magnification=abs(1./(1-poisson+DOUBLE(wl_dax_dx[xi:xf,yi:yf]*wl_day_dy[xi:xf,yi:yf]-wl_dax_dy[xi:xf,yi:yf]*wl_day_dx[xi:xf,yi:yf])*scaling_factor ))

				;find pixels where the lensed image ends up
				source_dist=sqrt((WLX[xi:xf,yi:yf]-wl_alpha_x[xi:xf,yi:yf]*scaling_factor-source_x)^2.+(WLY[xi:xf,yi:yf]-wl_alpha_y[xi:xf,yi:yf]*scaling_factor-source_y)^2.)
				indXY = WHERE(source_dist LT searchlensed)

					IF indXY[0] NE -1 THEN BEGIN
					;i.e. IF we find an image we save its position (indX,indY) and magnification mu
						indXY = ij(DOUBLE(indXY), N_ELEMENTS(poisson[*,0]))
						indX  = indXY[*,0]
						indY  = indXY[*,1]
						mu = magnification(indX,indY)
						indX  = indXY[*,0] + xi
						indY  = indXY[*,1] + yi
					ENDIF ELSE BEGIN
					;i.e. if we are here it means there is no image (probably it lies outside of the map)
					;we assign a non-sense position and magnification = 0.
						indX  = -999999.0
						indY  = -999999.0
						mu    = 0.0
					ENDELSE

					;we estimate the MEAN position (indX,indY) and magnification mu, in case the image was found in more than 1 pixel.
					mindX = mean(indX)
					mindY = mean(indY)
					mmu   = mean(mu)
				ENDELSE

;we have just finished the section that deals with sources that lie within the SL region ------------------------------------------------------
;the following section deals with sources that lie in the WL regime (very similar to the SL sections, but easier) -----------------------------

			ENDIF ELSE BEGIN
				print, 'WL regime'          ;just to know what the code is doing
				source_x = xx[ii]			;this is just the X position of the source in pixel units
				source_y = yy[ii]			;this is just the Y position of the source in pixel units

				;The next 4 IFs are to define the position of the WLmask to search for the source image.
				;The image will tend to appear outwards the cluter center (at least in the WL regime).
				;The mask is therefore not centered on the source position, but shifted depending on the quadrant where the source is.
				IF dfxx[ii] LE 0.0 AND dfyy[ii] LE 0.0 THEN BEGIN
					xi = source_x-WLmask
					xf = source_x+WLmask/2.
					yi = source_y-WLmask
					yf = source_y+WLmask/2.
				ENDIF

				IF dfxx[ii] LT 0.0 AND dfyy[ii] GT 0.0 THEN BEGIN
					xi = source_x-WLmask
					xf = source_x+WLmask/2.
					yi = source_y-WLmask/2.
					yf = source_y+WLmask
				ENDIF

				IF dfxx[ii] GT 0.0 AND dfyy[ii] GT 0.0 THEN BEGIN
					xi = source_x-WLmask/2.
					xf = source_x+WLmask
					yi = source_y-WLmask/2.
					yf = source_y+WLmask
				ENDIF

				IF dfxx[ii] GT 0.0 AND dfyy[ii] LT 0.0 THEN BEGIN
					xi = source_x-WLmask/2.
					xf = source_x+WLmask
					yi = source_y-WLmask
					yf = source_y+WLmask/2.
				ENDIF

				;These 4 IFs are just to make sure the WLmask do not fall outside the map
				IF xi LE 0.0 THEN xi = 0.0
				IF yi LE 0.0 THEN yi = 0.0
				IF xf GE N_ELEMENTS(wl_dax_dx[*,0])-1. THEN xf = N_ELEMENTS(wl_dax_dx[*,0])-1.
				IF yf GE N_ELEMENTS(wl_dax_dx[0,*])-1. THEN yf = N_ELEMENTS(wl_dax_dx[0,*])-1.

				;estimate magnifications (these come from lensing equations)
				poisson = DOUBLE(wl_dax_dx[xi:xf,yi:yf] + wl_day_dy[xi:xf,yi:yf])*scaling_factor
				magnification=abs(1./(1-poisson+DOUBLE(wl_dax_dx[xi:xf,yi:yf]*wl_day_dy[xi:xf,yi:yf]-wl_dax_dy[xi:xf,yi:yf]*wl_day_dx[xi:xf,yi:yf])*scaling_factor ))
				;find pixels where the lensed image ends up
				source_dist=sqrt((WLX[xi:xf,yi:yf]-wl_alpha_x[xi:xf,yi:yf]*scaling_factor-source_x)^2.+(WLY[xi:xf,yi:yf]-wl_alpha_y[xi:xf,yi:yf]*scaling_factor-source_y)^2.)
				indXY = WHERE(source_dist LT searchlensed)

				IF indXY[0] NE -1 THEN BEGIN
				;i.e. IF we find an image we save its position (indX,indY) and magnification mu
					indXY = ij(DOUBLE(indXY), N_ELEMENTS(poisson[*,0]))
					indX  = indXY[*,0]
					indY  = indXY[*,1]
					mu = magnification(indX,indY)
					indX  = indXY[*,0] + xi
					indY  = indXY[*,1] + yi
				ENDIF ELSE BEGIN
				;i.e. if we are here it means there is no image (probably it lies outside of the map)
				;we assign a non-sense position and magnification = 0.
					indX  = -999999.0
					indY  = -999999.0
					mu    = 0.0
				ENDELSE

				;we estimate the MEAN position (indX,indY) and magnification mu, in case the image was found in more than 1 pixel.
				mindX = mean(indX)
				mindY = mean(indY)
				mmu   = mean(mu)

			ENDELSE

;we have finished the WL regime section -------------------------------------------------------------------------------
		ENDELSE
;we have finished the lensing analysis section ------------------------------------------------------------------------

;finally!!!! here we construct the output arrays ------------------------------------------------------------

		IF ii EQ 0 THEN BEGIN

		x_out = mindX
		y_out = mindY
		mu_out = mmu

		id_in = dindgen(N_ELEMENTS(mindX))
		id_in[*] = dfid[ii]
		x_in = fltarr(N_ELEMENTS(mindX))
		x_in[*] = dfxx[ii]
		y_in = fltarr(N_ELEMENTS(mindX))
		y_in[*] = dfyy[ii]
		z_in = fltarr(N_ELEMENTS(mindX))
		z_in[*] = dfredz[ii]
		f_in = fltarr(N_ELEMENTS(mindX))
		f_in[*] = dfflux[ii]

		ENDIF ELSE BEGIN

		x_temp = dblarr(N_ELEMENTS(mindX) + N_ELEMENTS(x_out))
		x_temp[0:N_ELEMENTS(x_out)-1.] = x_out
		x_temp[N_ELEMENTS(x_out):*] = mindX
		x_out = x_temp

		y_temp = dblarr(N_ELEMENTS(mindY) + N_ELEMENTS(y_out))
		y_temp[0:N_ELEMENTS(y_out)-1.] = y_out
		y_temp[N_ELEMENTS(y_out):*] = mindY
		y_out = y_temp

		mu_temp = fltarr(N_ELEMENTS(mmu) + N_ELEMENTS(mu_out))
		mu_temp[0:N_ELEMENTS(mu_out)-1.] = mu_out
		mu_temp[N_ELEMENTS(mu_out):*] = mmu
		mu_out = mu_temp

		id_temp = dindgen(N_ELEMENTS(mindX) + N_ELEMENTS(id_in))
		id_temp[0:N_ELEMENTS(id_in)-1.] = id_in
		id_temp[N_ELEMENTS(id_in):*] = dfid[ii]
		id_in = id_temp

		xi_temp = dblarr(N_ELEMENTS(mindX) + N_ELEMENTS(x_in))
		xi_temp[0:N_ELEMENTS(x_in)-1.] = x_in
		xi_temp[N_ELEMENTS(x_in):*] = dfxx[ii]
		x_in = xi_temp

		yi_temp = dblarr(N_ELEMENTS(mindX) + N_ELEMENTS(y_in))
		yi_temp[0:N_ELEMENTS(y_in)-1.] = y_in
		yi_temp[N_ELEMENTS(y_in):*] = dfyy[ii]
		y_in = yi_temp

		zi_temp = dblarr(N_ELEMENTS(mindX) + N_ELEMENTS(z_in))
		zi_temp[0:N_ELEMENTS(z_in)-1.] = z_in
		zi_temp[N_ELEMENTS(z_in):*] = dfredz[ii]
		z_in = zi_temp

		fi_temp = dblarr(N_ELEMENTS(mindX) + N_ELEMENTS(f_in))
		fi_temp[0:N_ELEMENTS(f_in)-1.] = f_in
		fi_temp[N_ELEMENTS(f_in):*] = dfflux[ii]
		f_in = fi_temp

		ENDELSE

		;just to check how long it took to lens all the catalogue
		tf = systime(/seconds)
		print, 'time 2 = ', (tf-ti)/60.
	ENDFOR

;we have finished lensing the whole catalogue and we are ready to output resulte!
RETURN, TRANSPOSE([[id_in],[x_in],[y_in],[z_in],[f_in],[x_out],[y_out],[mu_out]])

END
;This is the end
;Beautiful friend
;This is the end
;My only friend, the end
