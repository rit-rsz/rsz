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
import config
from scipy.ndimage.filters import gaussian_filter
from scipy.stats import norm, chi2
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms

def confidence_ellipse(x, y, ax, n_std=1.0, facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    Returns
    -------
    matplotlib.patches.Ellipse

    Other parameters
    ----------------
    kwargs : `~matplotlib.patches.Patch` properties
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)

    return ax.add_patch(ellipse)

def chi_square_test(data,model,sigma):
    total_chi = []
    for i in range(3):
        total_chi.append((data[i] - model[i])**2 / (2*sigma[i]**2))
    # print('total_chi:',np.sum(total_chi),total_chi)
    return np.sum(total_chi)

nsim = 1
names = 'rxj1347'

# if sys.argv[1] == 'real' : # run the real maps
#     c(names, isim = None, saveplot=1,maketf=0,sgen=None,verbose=1,resolution='nr',superplot=0,testflag=0)
# else : # run the simulated maps
#     c(names, isim = sys.argv[1] ,saveplot=1,maketf=0,sgen=3,verbose=1,resolution='nr',superplot=0,testflag=1)
#
# exit()
#yt_grid and sz_grid should be indexed the same
# if os.path.isfile(config.HOME + 'outputs/%s_sz_grid.npy' %(names)):
#     sz_grid = np.load(config.HOME + 'outputs/%s_sz_grid.npy' %(names),allow_pickle=True)
#     input_yt = np.load(config.HOME + 'outputs/%s_input_grid.npy' %(names),allow_pickle=True)
#     ys = np.load(config.HOME + 'outputs/%s_y.npy' %(names))
#     ts = np.load(config.HOME + 'outputs/%s_t.npy' %(names))
# else :
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
like_li = []
sz_grid_0 = [x.get(0) for x in sz_grid]
sz_grid_1 = [x.get(1) for x in sz_grid]
sz_grid_2 = [x.get(2) for x in sz_grid]
rand_sz_grid_0 = np.random.normal(loc=sz_grid_0[int(len(sz_grid_0)/2.0)], scale=np.std(sz_grid_0), size=len(sz_grid_0))
rand_sz_grid_1 = np.random.normal(loc=sz_grid_1[int(len(sz_grid_1)/2.0)], scale=np.std(sz_grid_1), size=len(sz_grid_1))
rand_sz_grid_2 = np.random.normal(loc=sz_grid_2[int(len(sz_grid_2)/2.0)], scale=np.std(sz_grid_2), size=len(sz_grid_2))
# make chi square test for each final sz amplitude
# use sz_fit_real for real DIs
for i in range(len(sz_grid)):
    # chi_stat = chi_square_test([sz_fit_real[0][0],sz_fit_real[1][0],sz_fit_real[2][0]],
    #                             [sz_grid[i].get(0),sz_grid[i].get(1),sz_grid[i].get(2)],
    #                             [bias[0],bias[1],bias[2]])
    chi_stat = chi_square_test([sz_grid[i].get(0),sz_grid[i].get(1),sz_grid[i].get(2)],
                                [rand_sz_grid_0[i],rand_sz_grid_1[i],rand_sz_grid_2[i]],
                                [np.std(sz_grid_0),np.std(sz_grid_1),np.std(sz_grid_2)])

    like_li.append(np.exp(chi_stat/-2.0))

# compute likelihood function
YS,TS = np.meshgrid(ys,ts)
DI = np.array(like_li).reshape(len(ys),len(ts))
plt.pcolormesh(YS,TS,DI)
plt.colorbar()

# compute the 68.3% confidence contours
q = 2 * norm.cdf(1) - 1
r2 = chi2.ppf(q, 2)
print(q,r2)
cs = plt.contour(YS,TS,DI,levels=[q],colors=('k',),linestyles=('-',),linewidths=(2,))
plt.clabel(cs,fmt='%.3f',colors='k',fontsize=14)
plt.xlabel('Compton Y')
plt.ylabel('Temperature [K]')
plt.savefig('likelihood.png')
plt.clf()
