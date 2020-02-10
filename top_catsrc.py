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
import sys
sys.path.append('utilities')
sys.path.append('source_handling')
sys.path.append('reduc')
sys.path.append('sz')
from catsrc import Catsrc as c
from clus_szgf import *
from clus_sim_hist import *
from scipy.stats import chisquare, chi2
import config
from scipy.ndimage.filters import gaussian_filter

nsim = 200
length = 100
sim = [nsim + i for i in range(length)]
# names = 'a0370'
names = 'rxj1347'

for i in range(1):
    # c(names, saveplot=1,maketf=0,sgen=None,verbose=1,resolution='nr',superplot=0,testflag=1)
    # for j in range(100):
        #c(names, saveplot=1,maketf=0,sgen=2,nsim=sim[j],verbose=1,resolution='nr',superplot=0,testflag=0)
    c(names, saveplot=1,maketf=0,sgen=2,nsim=100,verbose=1,resolution='nr',superplot=0,testflag=1)

#yt_grid and sz_grid should be indexed the same
if os.path.isfile(config.HOME + 'outputs/%s_sz_grid.npy' %(names)):
    sz_grid = np.load(config.HOME + 'outputs/%s_sz_grid.npy' %(names),allow_pickle=True)
    input_yt = np.load(config.HOME + 'outputs/%s_input_grid.npy' %(names),allow_pickle=True)
    ys = np.load(config.HOME + 'outputs/%s_y.npy' %(names))
    ts = np.load(config.HOME + 'outputs/%s_t.npy' %(names))
else :
    sz_grid, input_yt, ys, ts = clus_szgf() # input_yt should have one DI for each band

avg_dI = clus_sim_hist(nsim,names)
print('avg dI : ',avg_dI)
print('input dI :',input_yt[0])

# calculate sz bias from pipeline
bias = [0]*3
for i in range(3):
    if avg_dI[i] < 0 :
        bias[i] = input_yt[0].get(i) + abs(avg_dI[i])
    else :
        bias[i] = avg_dI[i] - input_yt[0].get(i)
print('bias in dI : ',bias[0],bias[1],bias[2])
# what is the critical p value ?
crit_p = chi2.ppf(q = 0.683, df = 2)
print('critical value 68.3% confidence :',crit_p)

# grab real sz fit params
sz_fit_real = np.load(config.HOME + 'outputs/%s_fit.npy'%(names))
print(sz_fit_real)

# retreive and separate all dI for each band and y/t pair
# subtract bias for each band ... not sure that we have to do it by band
# maybe only if each bias is very different... we shall see
chisq_psw = []
chisq_pmw = []
chisq_plw = []
like_psw = []
like_pmw = []
like_plw = []
# make chi square test for each final sz amplitude
# use sz_fit_real for real DIs
for i in range(len(sz_grid)):
    stat1,pval1 = chisquare(sz_grid[i].get(0),f_exp=sz_fit_real[0][0] - bias[0],ddof=0)
    stat2,pval2 = chisquare(sz_grid[i].get(1),f_exp=sz_fit_real[1][0] - bias[1],ddof=0)
    stat3,pval3 = chisquare(sz_grid[i].get(2),f_exp=sz_fit_real[2][0] - bias[2],ddof=0)
    chisq_psw.append(stat1)
    like_psw.append(np.exp(stat1/-2.0))
    chisq_pmw.append(stat2)
    like_pmw.append(np.exp(stat2/-2.0))
    chisq_plw.append(stat3)
    like_plw.append(np.exp(stat3/-2.0))

# compute likelihood function
YS,TS = np.meshgrid(ys,ts)
DI = np.array(like_plw).reshape(len(ys),len(ts))
plt.pcolormesh(YS,TS,DI)
plt.colorbar()
cs = plt.contour(YS,TS,DI,levels=[1.172,1.173],colors=('k',),linestyles=('-',),linewidths=(2,))
plt.clabel(cs,fmt='%.3f',colors='k',fontsize=14)
plt.xlabel('Compton Y')
plt.ylabel('Temperature [K]')
plt.savefig('likelihood.png')
plt.clf()
