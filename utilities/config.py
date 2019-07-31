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
CLUSHOME = '/home/mercado/bitten/SPIRE/'
HOME = '/home/butler/rsz/'
FITSOUT = 'fits_images/'

calfac  = (pi/180.0)**2 * (1/3600.0)**2 * (pi / (4.0 * log(2.0))) * (1e6)
JY2MJy = 1e6 # Janskys to Mega Janskys

yin_coeff = [2.50,1.91,2.26,3.99,1.36,2.42,1.59,1.90,3.99]
yin = [x*1e-4 for x in yin_coeff]

tin = [7.2,10.1,7.7,9.8,4.5,8.6,7.8,5.5,10.9]