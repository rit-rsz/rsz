################################################################################
# NAME : save_data.py
# DATE STARTED : June 24, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE : idk lol
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import numpy as np

def save_data(maps, yin=0, tin=0, simflag=0,verbose=1):
    ncol = len(maps)
    for i in range(ncol):
        thismap = maps[i]
        if thismap['band'] in thismap['file']# need to figure out what to put here :)
            ststart = thismap['file'].find(thismap['band'])
            stend = thismap['file'].find('.')
            fileext = thismap['file'][ststart:stend] + '.npy'

        if simflag == 0:
            filename = config.CLUSDATA + 'sz/' + thismap['name'] + '_' + fileext
            print(filename)

        else:
            thisy = str(yin / (1*10**-4))
            thist = str(tin)
            filename = config.CLUSDATA + 'sz/sim/' + thismap['name'] + '/' + thismap['name'] + '_' + thisy + '_' thist + '_' + fileext

        #not sure what the point of the if statements are if its just going to set the filename to thismap['file'] anyway...
        thismap['file'] = filename
        np.save(filename, thismap) #i have a feeling thismap needs to be an array.
        #call to this into a .sav file currently looking for a python equivalent.
