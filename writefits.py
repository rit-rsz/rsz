################################################################################
# NAME : writefits.py
# DATE STARTED : June 24, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE : to create a fits file with the given information
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS : filename - the name of the file
#          header_dict - an optional argument that contains the header information
#          data - an optional argument that contains the data.
#
#
# OUTPUTS : None - but it creates a file
# REVISION HISTORY :
################################################################################


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
