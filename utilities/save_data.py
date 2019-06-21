################################################################################
# NAME : save_data.py
# DATE STARTED : June 21, 2019
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

def save_data(maps, yin=0, tin=0, simflag=0,verbose=1):
    ncol = len(maps)

    for i in range(ncol):
        thismap = maps[i]
        if thismap['band'] in thismap['file']# need to figure out what to put here :)
        
