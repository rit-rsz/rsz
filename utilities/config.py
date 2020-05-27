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

CLUSDATA = '/data/butler/SPIRE/'
# CLUSSBOX = '/data/butler/SPIRE/sandbox/'
CLUSHOME = '/home/vlb9398/bitten/SPIRE/'
HOME = '/home/butler/rsz/'
#SIM = '/home/vlb9398/rsz/new_bethermin/'
#CLUSSIMS = '/home/vlb9398/data/SPIRE/bethermin_sims/'
ROOT = '/home/vlb9398/'
OUTPUT = '/home/vlb9398/rsz/outputs/'
CLUSSIMS = '/data/butler/SPIRE/new_bethermin_sims/'
# CLUSSIDES = '/data/butler/SPIRE/sides_sims/'
SIMBOX = '/home/butler/rsz/new_bethermin/lens_model/'

calfac  = (pi/180.0)**2 * (1/3600.0)**2 * (pi / (4.0 * log(2.0))) * (1e6)
JY2MJy = 1e6 # Janskys to Mega Janskys

yin_coeff = [4.56,8.5,8.7,4.45,7.43,8.96,0.93,7.42,5.63,3.53,3.72,14.08]
yin = [x*1e-4 for x in yin_coeff]

tin = [7.2,10.1,7.65,6.7,9.81,9.16,4.5,8.62,7.8,7.2,5.5,10.88]

''' These were the ones used in the IDL pipeline '''
# yin_coeff = [2.50,1.91,2.26,3.99,1.36,2.42,1.59,1.90,14.08]
# yin = [x*1e-4 for x in yin_coeff]
#
# tin = [7.2,10.1,7.7,9.8,4.5,8.6,7.8,5.5,10.88]
