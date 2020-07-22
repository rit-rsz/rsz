################################################################################
# Name : get_iris_map
# Purpose : this script contains a function that is meant to retrieve
# information from the IRIS maps.
#Author : Benjamin Vaughan
#Start Date : October 4, 2019
#Additional Info
#
################################################################################
import os
from astropy.io import fits

def get_iris(ii, dir='!ISRISDATA', band=4, hcon=0, verbose=1):
    '''
    Purpose : returns the correlated IRIS map for an input ISSA number
    Inputs  : ii - the ISSA map number
              verbose - verbosity flag
              band - The IRAS band
                - 1 : 12 micron
                - 2 : 25 micron
                - 3 : 60 micron
                - 4 : 100 micron (default)
              hcon - hcon number (default is zero -> co-added map)
              dir - directory where the IRIS data is stored (default is !IRISDATA in idl)
    Outputs : map - the ISSA map corresponding to number ii
    '''
    iras_number = str(int(ii))
    if ii < 10:
        iras_number = '0' + iras_number
    elif ii < 100:
        iras_number = '0' + iras_number

    files = []
    for x in os.listdir(dir):
        if iras_number in x and str(band) in x and str(hcon) in x:
            files.append(dir + '/' + x)
    if len(files) > 0:
        hdul = fits.open(files[0])

        header = hdul[0].header
        #force lonpole to be 180
        try:
            del(header['LONPOLE'])
        except KeyError:
            pass
        header.append(('LONPOLE',180))
        #the third dimension on the IRIS fits files is only 1 element and Astropy's
        #built in pixel <-> skycoord transformation functions will try to do a 3D
        #transformation rather than 2D with this extra dimension in the header.
        if header['CDELT3'] == 0 and header['NAXIS3'] == 1:
            del(header['CDELT3'])
            del(header['NAXIS3'])
            del(header['CRPIX3'])
            del(header['CRVAL3'])
            del(header['CTYPE3'])

            header['NAXIS'] = 2

        #check data for bad values
        for i in range(hdul[0].data.shape[0]):
            for j in range(hdul[0].data.shape[1]):
                #we have to add this extra zero because the IRIS maps are
                #3 Dimensional with a z axis that is 1 element thick.
                if hdul[0].data[i,j,0] < -5 or hdul[0].data[i,j,0] == 0:
                    hdul[0].data[i,j] = -32768
        hdu = fits.PrimaryHDU(hdul[0].data, header)
        map = fits.HDUList([hdu])
        if verbose:
            print('Read data file %s' % files[0])
        return map
    else:
        print('Could not find any files matching that description')
        return None
