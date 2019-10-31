################################################################################
# NAME : test.py
# DATE STARTED : October 25, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE : This code is just for a first test run of actually running through
# the whole data analysis pipeline and returning graphs of the results.
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
from catsrc import Catsrc as c



nsim = 200
length = 100
sim = [nsim + i for i in range(length)]

x = c('a0370', saveplot=1, nsim=0, verbose=0, superplot=0, cattype='PSW')
print('successful')
for i in range(100):
    x = c('a0370', saveplot=1, simmap=2, nsim=sim[i], verbose=0, superplot=0, cattype='PSW')
