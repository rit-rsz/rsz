from astropy.io import fits

def writefits(dict):
    header = fits.Header(dict, filename, data=none)
    hdu = fits.PrimaryHDU(header=header, data=data)
    hdu.writeto(filename)
