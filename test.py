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
t = [7.2, 10.1, 10.9, 5.5]
y = [2.5, 1.91, 3.99, 1.9]
yin = [yi * 10**-4 for yi in y]
names = ['a0370', 'a1689', 'rxj1347', 'rxj0152']

#ben
for i in range(0, 1, 1):
    x = c(names[i], saveplot=1, nsim=0, verbose=0, superplot=0, saveplot=1, cattype='PSW', yin=yin[i], tin=t[i])
    for j in range(100):
        x = c(names[i], saveplot=1, simmap=2, nsim=sim[j], verbose=0, superplot=0, saveplot=1, cattype='PSW', yin=yin[i], tin=t[i])

# #dale
# for i in range(2, 3, 1):
#     x = c(names[i], saveplot=1, nsim=0, verbose=0, superplot=0, saveplot=1, cattype='PSW', yin=yin[i], tin=t[i])
#     for j in range(100):
#         x = c(names[i], saveplot=1, simmap=2, nsim=sim[j], verbose=0, superplot=0, saveplot=1, cattype='PSW', yin=yin[i], tin=t[i])
