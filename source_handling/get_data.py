
################################################################################
# NAME : get_data.py
# DATE STARTED : June 11, 2019
# AUTHORS : Dale Mercado & Benjamin Vaughan
# PURPOSE : This is where the data is retevied and can be
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import scipy.io
import numpy as np
# from config import * #(this line will give me access to all directory variables)
import matplotlib.pyplot as plt
import math
from astropy.io import fits
import os
import sys
import pyfits
import config
sys.path.append('utilities')
from get_spire_beam import *

def get_data(clusname, verbose = 1, resolution = 'nr', version = '1', manpath=0,
             bolocam=None, manidentifier=None):
    # place holders for the if statements to work until I add them to the input for get data
    # This will happen when the script is fully functional
    errmsg = False
    if not bolocam:
        cols = ['PSW','PMW','PLW']
        bolocam = 0
    else:
        cols = ['PSW','PMW','PLW','BOLOCAM']
        bolocam = 1

#   If there is no manual path set
    if manpath == 0:
        hermfiles = []
        hermdir = config.CLUSDATA + 'hermes_clusters/'
        hermcount = 0
        for x in os.listdir(hermdir):
            if x.startswith(clusname):
                if resolution in x and version in x:
                    hermfiles.append(hermdir + x) #this SHOULD work
                    hermcount += 1
        hlsdir = config.CLUSDATA + 'hls_clusters/'
        hlsfiles = []
        hlscount = 0
        for x in os.listdir(hlsdir):
            if x.startswith(clusname):
                if x.startswith(clusname):
                    hlsfiles.append(hlsdir + x)
                    hlscount += 1
        snapdir = config.CLUSDATA + 'snapshot_clusters'
        snapfiles = []
        snapcount = 0
        for x in os.listdir(snapdir):
            if x.startswith(clusname):
                snapfiles.append(snapdir + x)
                snapcount += 1


        # hermfile = [self.dir + x for x in x.startswith('/data/mercado/SPIRE/' + 'hermes_cluster/' + clusname)\
        #             + x.endswith('_' + resolution + '_' + version + '.fits')]
        # files = [self.dir + x for x in os.listdir(self.dir) if x.startswith("omnilog")]
        # This is still broken and I am usure how to move forward on this
        # print(hermfile)
        # hermfile = str('/data/mercado/SPIRE/' + 'hermes_clusters/' + clusname + '*_' +\
        #             resolution + '_' + version + '.fits')
        # This is an example of the file name that we are trying to call
        # rxj1347_PLW_fr_1.fits
        # hermfiles = []
#        hermcount = 0
#        hlscount = 0
#        snapcount = 0
#        for i in range(len(cols)):
            # Still need to figureout how to append these fits files to an array of some sort
#            hermfile = (CLUSDATA + 'hermes_clusters/' + clusname + '_' +\
#                        cols[i] + '_' + resolution + '_' + version + '.fits')
#            hermfiles = (fits.open(hermfile))

#           Count the number of hermfiles that are looked at
#            if len(hermfiles):
#                hermcount += 1
#            print(hermcount)

            # In the idl file it uses a * to grab all the bands here I had to use a for loop in order
            # to acheive the same results
            # This method provides difficulty for the hls and snapfiles below as they use other numbers in their names


#           Currently these two are broken:
#           I need to figureout how to use the * to call them all.
#           There are other variables not that are being called and I cant work around it like the otehr
#
#           This is example of all the hls files that are being called
#           a0068_BL_HP16_1.0_8.1.fits
#           a0068_PLW_12_8.2.fits
#           a0068_PMW_9_8.2.fits
#           a0068_PSW_6_8.2.fits
#           a0068_R_HP25_2.0_8.1.fits

            # hlsfile = (CLUSDATA + 'hls_clusters/' + clusname + '*')
            # hlsfiles = fits.open(hlsfile)
            # print(hlsfiles.info())

#           This is an example of all of the snapshot files that are getting called,
#           only 3 in the dicrectory
#           cl1226_PLW_12_8.2.fits
#           cl1226_PMW_9_8.2.fits
#           cl1226_PSW_6_8.2.fits
            # snapfile = (CLUSDATA + 'snapshot_clusters/' + clusname +'*.fits')
            # snapfiles = fits.open(snapfile)

#       For now I will just create and empty array for these ones and come back when I havea better idea of how to proceed


        if snapcount ^ hermcount ^ hlscount != 3 :
            errmsg = ('Problem finding exclusive files, file dump is:', \
                        snapfiles, hermfiles, hlsfiles)
            if verbose:
                print(errmsg)
            return None, errmsg
        if hermcount == 3:
            files = hermfiles
            nfiles = hermcount

        if hlscount == 3:
            files = hlsfiles
            nfiles = hlscount

        if snapcount == 3:
            files = snapfiles
            nfiles = snapcount
        print('help') #please help me

#   Manual Path option
    else:
        mancount = 0
        manfiles = []
        file_flag = False
        for x in os.listdir(manpath):
            if clusname in x and manidentifier in x:
                manfiles.append(manpath + x)

        if len(manfiles) != 3:
            errmsg = 'Cannot find files in %s' %(manpath)
            if verbose:
                print(errmsg)
            return None, errmsg
        nfiles = len(manfiles)
        files = manfiles


    if bolocam:
        nfiles = nfiles +1


#   This is maps as the ptr array<NullPointer><NullPointer><NullPointer>
#   nfiles = 3 I assume that accounts for the three null pointers being populated
    maps = []

#   Need to tweek the syntax of this for loop
    for ifile in range(nfiles):
        count = []
        counter = 0
        if ifile < 3:
            if cols[ifile] in files[ifile]:
                counter += 1
            else:
                errmsg = 'Problem finding ' + cols[ifile] + ' file.'
                if verbose:
                    print(errmsg)
            maps.append(read_file(files[counter], cols[ifile], clusname, verbose=verbose))
        else:
                maps[ifile] = np.empty(clus_read_bolocam(clusname,verbose=verbose)) #args need to be filled in bolocam one
    return maps, errmsg

#               From get_clus params it's very close to the same
# param_len = len(param_data)
# count = 0
# for i in range(param_len):
#     place = str(param_name[i])
#
#     if place == clusname_in :
#         whpl = i
#         count += 1

#           makes maps a pointer again, we're gona try using a dictionary
#           This is also where readfile is used for the first time
        #     maps[ifile] = read_file(files[whpl],cols[ifile],\
        #                                   clusname,VERBOSE=verbose)
        # else:
        #     maps[ifile] = read_file(files[whpl],cols[ifile],\
        #                                   clusname,VERBOSE=vercolbose)



def read_file(file,col,clusname,verbose=0):
#   The way this works in IDL is a pointer is populated with this kind of information in this case
#   we should be using a dictionary as this is one of the important ones called maps
#   This should ulitmatly be defined in as a self.calfac for this calibration
    calfac = (np.pi/180) * (1/3600)**2 * (np.pi / (4 * np.log(2))) * (1e6)
#   This will ultimatly be in the list of constants
    # The rest of the scrpit involves idl_libs stuff that
    # will get grabbed from astropy
    hdul = fits.open(file)
    map = hdul[1]
    err = hdul[2]
    exp = hdul[3]
    flag = hdul[4]

    if 'CDELT1' in map.header.keys():
        pixsize = 3600 * mean([abs(map[0].header['CDELT1'],abs(map[0].header['CDELT2']))])
        map[0].header['cd_1'] = map[0].header['CDELT1']
        map[0].header['cd_2'] = 0
        map[0].header['cd_1'] = 0
        map[0].header['cd_2'] = map[0].header['CDELT2']
        err[0].header['cd_1'] = err[0].header['CDELT1']
        err[0].header['cd_2'] = 0
        err[0].header['cd_1'] = 0
        err[0].header['cd_2'] = err[0].header['CDELT2']

    else:
        pixsize = 3600 * \
                  mean([abs(map.header['CD1_1']+map.header['CD2_1']), \
                        abs(map.header['CD2_1'] + map.header['CD2_2'])])

    psf = get_spire_beam(pixsize, band=col, factor=1)
    widtha = get_spire_beam_fwhm(col)
    width = (widtha / math.sqrt(8 * math.log(2)) * pixsize)
#   We wouldnt be able to put this one in calfac since it is determined by the source called
    calfac = 1 / (calfac * (get_spire_beam_fwhm(col))**2)
#   This should be defined in the catsrsc file
    JY2MJy = 1e6

#   Something called EXTAST?????    astr = map.header

#   Gets header information from a fits image. Astropy should be able to do this
    astr = {}
    try:
        cd11 = map.header['CD1_1']
        cd12 = map.header['CD1_2']
        cd21 = map.header['CD2_1']
        cd22 = map.header['CD2_2']
        pv1_1 = map.header['PV1_0']
        pv1_2 = map.header['PV1_1']
        pv1_3 = map.header['PV1_2']
        pv1_4 = map.header['PV1_3']
        pv1_5 = map.header['PV1_4']
    except KeyError:
        pass # i don't like the way this is coded probably have to change it later
    for keys in map.header.keys():
        if 'NAXIS' in keys:
            astr.update({keys : map.shape})
        if 'CD1_1' in keys:
            x =  np.array([[cd11, cd12], [cd21, cd22]])
            astr.update({'CD' : x})
        if 'CDELT' in keys:
            astr.update({keys : map.header[keys]})
        if 'CRPIX1' in keys:
            x = np.array([map.header['CRPIX1'], map.header['CRPIX2']])
            astr.update({'CRPIX' : x})
        if 'CTYPE1' in keys:
            x = np.array([map.header['CTYPE1'], map.header['CTYPE2']])
            astr.update({'CTYPE' : x})
        if 'CRVAL1' in keys:
            x = np.array([map.header['CRVAL1'], map.header['CRVAL2']])
            astr.update({'CRVAL' : x})
        if 'LONGPOLE' in keys:
            astr.update({keys : map.header[keys]})
        if 'LATPOLE' in keys:
            astr.update({keys : map.header[keys]})
        if 'PV1_0' in keys:
            x = np.array([pv1_1, pv1_2, pv1_3, pv1_4, pv1_5])
            astr.update({keys : x})
        if True:
            pass

    head = map.header
    herr = err.header

#   Not sure if this is the correct syntax for astr naxis
    srcrm = np.empty(astr['NAXIS'])
    xclean = np.empty(astr['NAXIS'])
    mask = np.empty(astr['NAXIS'])

#   generate default mask map
#   searches though the data and if the log = 0 the value is filled with NaN
#   it marks the place where if map is finite and = 0 with countnan
    # whnan = np.where(np.isfinite(map.data) == False)
    # if len(whnan) > 0:
    #     pass #i really have no idea what to put here i think i don't fundamentally understand what is happening in the idl script


#    This is how the idl to python site shows how to do these conditional indexing
	# (i,j) = a.nonzero()
    # (i,j) = where(a!=0)

    # SXPAR()  data[0].header['KEY'] 	Obtain value of header keyword
    # SXADDPAR()  data[0].header['KEY']='value' OR data[0].header.update('KEY','value', 'comment')

    maps = {'name':clusname, #check
          'file':file, #check
          'band':col, #check
          'signal':map, #check
          'srcrm':srcrm, #check
          'xclean':xclean, #check
          'error':err, #check
          'mask':mask, #check
          'flag':flag, #check
          'exp':exp, #check
          'shead':head, #nope
          'ehead':herr, #nope
          'astr':astr, #check
          'pixsize':pixsize, #check
          'psf':psf, #check
          'width':width, #check
          'widtha':widtha, #check
          'calfac':calfac, #check
          'JY2MJy':JY2MJy} #check

    return maps

def mean(data):
    length = len(data)
    sum = 0
    for i in range(length):
        sum += data[i]
    average = sum / length
    return average

if __name__ == '__main__':
    get_data('rxj1347',1)
