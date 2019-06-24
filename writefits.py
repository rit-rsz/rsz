from astropy.io import fits

def writefits(filename, header_dict=None, data=None):
    if header_dict:
        header = fits.Header(header_dict)
    if header and data:
        hdu = fits.PrimaryHDU(header=header, data=data)
    elif header:
        hdu = fits.PrimaryHDU(header=header)
    elif data:
        hdu = fits.PrimaryHDU(data=data)
    hdu.writeto(filename)
