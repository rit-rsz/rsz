import matplotlib.pyplot as plt
from genmap import genmap_gauss
import time
import numpy as np
from scipy.interpolate import UnivariateSpline
import math, sys

m250 = []
index = []
gm = genmap_gauss()
maps = gm.generate(0.25,verbose=True)
data = maps[-1]
fluxes = data['fluxdens']

for i in range(len(fluxes)):
    m250.append(fluxes[i][0] * 1e3)

# avg_flux = np.log10(avg_flux)
# for i in range(len(avg_flux)):
#     if avg_flux[i] < 0.0 :
#         index.append(i)

num_bin = math.ceil(max(m250)-min(m250)) - 1
# size_bin = 10 * np.linspace(0.0,np.log10(max(m250)),num_bin)
size_bin = np.linspace(1.0,max(m250),num_bin)
print('num_bin: ',num_bin)

# hist,bin = np.histogram(m250,bins= 10 * np.linspace(0.0,np.log10(max(m250)),num_bin))
hist,bin = np.histogram(m250,bins= np.linspace(1.0,max(m250),num_bin))
hist = [(float(x) / ((0.25/((180/np.pi)**2))*((size_bin[3] - size_bin[2]) / 1e3))) for x in hist]
print(np.diff(size_bin,n=1)[0:10])
for i in range(len(hist)):
	hist[i] = np.log10(hist[i] * ((size_bin[i]/1e3)**(2.5)))

# bins = bin[:-1] + (bin[1]-bin[0])/2
# f = UnivariateSpline(bins,hist)
# plt.scatter(bins,f(bins),color='orange')
# plt.plot(bins,f(bins),color='black')
# y_val = interp1d(size_bin,hist)
# plt.scatter(y_val(10.0))
plt.hist(bin[:-1],bin,weights=hist,histtype='step')
plt.gca().set_xscale("log")
# plt.gca().set_yscale("log")
# plt.xlim((0.1,1000))
# plt.ylim((0.1,1e10))
plt.xlabel('S [mJy]')
plt.ylabel('dN/dS [Jy$^1.5$/sr$^-1$]')
plt.title('250 microns')
plt.show()
