
############################################################################
#
# NAME : clus_format_bethermin.py
# DATE : October 3, 2019
# AUTHOR : Victoria Butler
# PURPOSE : takes the catalog output by Alex Conley's implementation of
#           the bethermin model and puts it into the format that lenstool wants.
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
# OUTPUTS :
# REVISION HISTORY :

############################################################################


def clus_format_bethermin(icol,catfile,mapfile,band,clusname,pixsize,\
                          FLUXCUT=fluxcut,ZZERO=zzero,RETCAT=retcat):

'''
  ; standard init stuff
  COMPILE_OPT IDL2, STRICTARRSUBS
  success = 0b
  IF ~(N_ELEMENTS(verbose)) THEN verbose=1 ELSE verbose = verbose

  IF ~N_ELEMENTS(fluxcut) THEN fluxcut=0 ELSE fluxcut=fluxcut
  IF ~N_ELEMENTS(zzero) THEN zzero=0 ELSE zzero = zzero

  ;; won't error check arguments because you'd better know what
  ;; you're doing to be playing with this program!

  msrc = 60000L - 1L

  cat = MRDFITS(catfile,1)
  nsrc = N_ELEMENTS(cat)

  mapfilename = STRING(mapfile,'_',band,'.fits')
  map = MRDFITS(mapfilename,1,header)

  refx = SXPAR(header,'CRPIX1')
  refy = SXPAR(header,'CRPIX2')

  ;; massage data into new arrays
  ;; this is going to break at pixel sizes different from 1"
  outx = pixsize*(cat.xpos - refx)
  outy = pixsize*(cat.ypos - refy)
  outz = FLOAT(ROUND(10. * cat.z)) / 10.
  outflux = cat[*].fluxes[icol]

  IF fluxcut GT 0 THEN BEGIN
     MESSAGE,'Cutting input catalog at flux ' + STRING(fluxcut),/INFORMATIONAL
     whpl = WHERE(outflux GT fluxcut,count)
     IF count GT 0 THEN BEGIN
        outx = outx[whpl]
        outy = outy[whpl]
        outz = outz[whpl]
        outflux = outflux[whpl]
     ENDIF
     nsrc = N_ELEMENTS(outflux)
  ENDIF
  IF zzero GT 0 THEN BEGIN
     whpl = WHERE(REFORM(outz) GT zzero[0],count,COMPLEMENT=loz)
     IF count GT 0 THEN BEGIN
        retcat = {x:outx[loz],y:outy[loz],z:outz[loz],f:outflux[loz]}
        outx = outx[whpl]
        outy = outy[whpl]
        outz = outz[whpl]
        outflux = outflux[whpl]
     ENDIF
     nsrc = N_ELEMENTS(outflux)
  ENDIF

  ;; sort according to brightness due to lenstool limitations
  sortf = REVERSE(SORT(outflux))
  outx = outx[sortf]
  outy = outy[sortf]
  outz = outz[sortf]
  outflux = outflux[sortf]

  IF msrc LT nsrc THEN BEGIN
     ;; truncacte to the msrc brightest sources
     outx = outx[0:msrc]
     outy = outy[0:msrc]
     outz = outz[0:msrc]
     outflux = outflux[0:msrc]

     nsrc = msrc
  ENDIF

;top

  ;; now sort according to z
  sortz = SORT(outz)
  outx = outx[sortz]
  outy = outy[sortz]
  outz = outz[sortz]
  outflux = outflux[sortz]

  outmag = -2.5 * ALOG10(outflux)

  lensfile = STRING(!clushome,'model/',clusname,'/',clusname,'_cat.cat')
  OPENW,1,lensfile
  mystring = '#REFERENCE 3 ' + $
             STRCOMPRESS(STRING($
             SXPAR(header,'CRVAL1'),FORMAT='(F10.6)'),/REMOVE_ALL) + ' ' + $
             STRCOMPRESS(STRING($
             SXPAR(header,'CRVAL2'),FORMAT='(F10.6)'),/REMOVE_ALL)
  PRINTF,1,mystring
  FOR icat=0,nsrc-1L DO BEGIN
     mystring = STRCOMPRESS(STRING(icat+1),/REMOVE_ALL) + '   ' + $
                STRCOMPRESS(STRING(outx[icat]),/REMOVE_ALL) + '   ' + $
                STRCOMPRESS(STRING(outy[icat]),/REMOVE_ALL) + '   ' + $
                '0.5   0.5   0.0   ' + $ ;2.5   ' + $
                STRCOMPRESS(STRING(outz[icat]),/REMOVE_ALL) + '   ' + $
                STRCOMPRESS(STRING(outmag[icat]),/REMOVE_ALL)
     PRINTF,1,mystring
  ENDFOR
  CLOSE,1


  RETURN

END
'''
