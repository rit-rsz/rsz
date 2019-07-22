################################################################################
# NAME : config.py
# DATE STARTED : May 30, 2019
# AUTHORS : Victoria Butler & Dale Mercado & Benjamin Vaughan
# PURPOSE : This file holds the directories and variables common to all scripts
#           run in this series
# EXPLANATION :
# CALLING SEQUENCE
# INPUTS :
# OUTPUTS :
# REVISION HISTORY :
################################################################################
from math import *

CLUSDATA = '/data/mercado/SPIRE/'
CLUSSBOX = '/data/mercado/SPIRE/sandbox/'
FITSOUT = 'fits_images/'

calfac  = (pi/180.0) * (1/3600.0)**2 * (pi / (4.0 * log(2.0))) * (1e6)
JY2MJy = 1e6 # Janskys to Mega Janskys
