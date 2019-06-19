
################################################################################
# NAME : get_data.py
# DATE STARTED : June 11, 2019
# AUTHORS : Dale Mercado
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
import pyfits
import config

def get_data(clusname, verbose = 1, resolution = 'nr', version = '1', manpath=0,
             bolocam=None):
    # place holders for the if statements to work until I add them to the input for get data
    # This will happen when the script is fully functional

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
            return success
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
        for x in os.listdir(manpath):
            pass # don't know what to put here because i don't know what the specifications for manpath are
        if len(files):
            mancount += 1
        if mancount == 3:
            errmsg = 'Cannot find 3 files fitting that path description!'
            #message errmsg
            print(errmsg)
            return success
        nfiles = mancount


    if bolocam:
        nfiles = nfiles +1


#   This is maps as the ptr array<NullPointer><NullPointer><NullPointer>
#   nfiles = 3 I assume that accounts for the three null pointers being populated
    maps = nfiles

#   Need to tweek the syntax of this for loop
    for ifile in range(nfiles):
        count = []
        counter = 0
        if ifile < 3:
            for i in range(len(nfiles)):
                counter += 1
                if cols[ifile] in files:
                    count.append(counter)
            if len(count) != 1:
                print('Problem finding', cols[ifile], 'file.')
            else:
                maps[ifile] = clus_read_file(args) #args need to be filled in
        else:
            maps[ifile] = clus_read_file(args) #args need to be filled in bolocam one
    return maps

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
        #                                   clusname,VERBOSE=verbose)



def read_file(file,col,clusname,VERBOSE=verbose):
#   The way this works in IDL is a pointer is populated with this kind of information in this case
#   we should be using a dictionary as this is one of the important ones called maps
#   This should ulitmatly be defined in as a self.calfac for this calibration
    calfac = (np.pi/180) * (1/3600)**2 * (np.pi / (4 * np.log(2))) * (1e6)
#   This will ultimatly be in the list of constants

    # The rest of the scrpit involves idl_libs stuff that
    # will get grabbed from astropy
    map = fits.open(file[0])
    err = fits.open(file[1])
    exp = fits.open(file[2])
    flag = fits.open(file[3])

    if map[0].header['CDELT1'] != 0:
        pixsize = 3600 * \
                  mean([abs(map[0].header['CDELT1'],abs(map[0].header['CDELT2']))])
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
                  mean([abs(map[0].header['CD1_1']+map[0].header['CD2_1']), \
                        abs(map[0].header['CD2_1'] + map[0].header['CD2_2'])])

    psf = get_spire_beam(col,pixsize)
    widtha = get_spire_beam_fwhm(col)
    width = widtha / sqrt(8 * log(2)) * pixsize)
#   We wouldnt be able to put this one in calfac since it is determined by the source called
    calfac = 1 / (calfac * (get_spire_beam_fwhm(col))**2)
#   This should be defined in the catsrsc file
    JY2MJy = 1e6

#   Something called EXTAST?????
#   Gets header information from a fits image. Astropy should be able to do this
    astr = map.header[]

#   Not sure if this is the correct syntax for astr naxis
    srcrm = np.zeros(astr.naxis)
    xclean = np.zeros(astr.naxis)
    mask = np.zeros(astr.naxis)

#   generate default mask map
#   searches though the data and if the log = 0 the value is filled with NaN
#   it marks the place where if map is finite and = 0 with countnan
    whnan = np.where(np.isfinite(map) == False)
    if len(whnan) > 0:
        mask[whnan] = 1

#    This is how the idl to python site shows how to do these conditional indexing
	# (i,j) = a.nonzero()
    # (i,j) = where(a!=0)

    # SXPAR()  data[0].header['KEY'] 	Obtain value of header keyword
    # SXADDPAR()  data[0].header['KEY']='value' OR data[0].header.update('KEY','value', 'comment')

    maps = {'name':clusname,\
          'file':file,\
          'band':col,\
          'signal':map,\
          'srcrm':srcrm,\
          'xclean':xclean,\
          'error':err,\
          'mask':mask,\
          'flag':flag,\
          'exp':exp,\
          'shead':head,\
          'ehead':herr,\
          'astr':astr,\
          'pixsize':pixsize,\
          'psf':psf,\
          'width':width,\
          'widtha':widtha,\
          'calfac':calfac,\
          'JY2MJy':JY2MJy}

    return maps

if __name__ == '__main__':
    get_data('rxj1347',1)
