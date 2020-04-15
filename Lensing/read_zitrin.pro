
FUNCTION READ_ZITRIN,npix,filename
OPENR,1,filename
Z_array = DBLARR(npix,npix)
READF,1,Z_array
CLOSE,1
RETURN, Z_array
END
