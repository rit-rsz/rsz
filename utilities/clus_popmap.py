
##########################################################################
# NAME : clus_popmap.py
# DATE STARTED : October 3, 2019
# AUTHORS : Victoria Butler
# PURPOSE : This program populates simulated maps with sources lensed by lenstool.
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :

# OUTPUTS :
# REVISION HISTORY :

##########################################################################


def clus_popmap(ltfile,copymap,MAP=outmap,LOZ=loz):

'''
  IF ~N_ELEMENTS(loz) THEN loz=0 ELSE loz=loz

  READCOL,ltfile,dummy,type,racent,deccent,NUMLINE=1,$
          FORMAT='A,A,F,F',/SILENT

  READCOL,ltfile,id,ra,dec,xs,ys,theta,z,mag,SKIPLINE=1;,/SILENT

  lozsize = SIZE(loz)
  IF lozsize[2] EQ 8 THEN BEGIN
     ra = [ra,loz.x]
     dec = [dec,loz.y]
  ENDIF

  ra = -ra / 3600. + REPLICATE(racent,N_ELEMENTS(ra))
  dec = dec /3600. + REPLICATE(deccent,N_ELEMENTS(dec))
  ;; this 4 shouldn't be here and should be removed if you rerun sims
  flux = 10.^(-mag / 2.5); / 4.
  IF lozsize[2] EQ 8 THEN BEGIN
     flux = [flux,loz.f]
  ENDIF

  AD2XY,ra,dec,copymap.astr,x,y

  psf = copymap.psf ;/ TOTAL(copymap.psf)

  outmap = FLTARR(copymap.astr.naxis)
  FOR isrc=0l,N_ELEMENTS(ra) - 1l DO BEGIN
     IF x[isrc] GE 0 and x[isrc] LT copymap.astr.naxis[0] AND $
        y[isrc] GE 0 AND y[isrc] LT copymap.astr.naxis[1] THEN $
           outmap[x[isrc],y[isrc]] = outmap[x[isrc],y[isrc]] + flux[isrc]
  ENDFOR

  outmap = CONVOLVE(outmap,psf)

  ;outmap = IMAGE_MODEL(x,y,flux,copymap.astr.naxis[0],copymap.astr.naxis[1],psf)

  outfile = STRING(!CLUSSBOX,'lensedmap_',copymap.name,'_',copymap.band,'.fits')
  WRITEFITS,outfile,outmap,copymap.shead

  RETURN

END
'''
