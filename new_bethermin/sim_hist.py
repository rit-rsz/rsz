import matplotlib.pyplot as plt
from genmap import genmap_gauss
import time
import numpy as np
from scipy.interpolate import UnivariateSpline
import math, sys

m250 = []
m350 = []
m500 = []
band = ['250 microns','350 microns','500 microns']
gm = genmap_gauss()
maps = gm.generate(0.25,verbose=True)
data = maps[-1]
fluxes = data['fluxdens']

for i in range(len(fluxes)):
    m250.append(fluxes[i][0] * 1e3)
    m350.append(fluxes[i][1] * 1e3)
    m500.append(fluxes[i][2] * 1e3)

big_data = [m250,m350,m500]
for i in range(len(big_data)):
    num_bin = math.ceil(max(big_data[i])-min(big_data[i])) - 1
    # size_bin = 10 * np.linspace(0.0,np.log10(max(m250)),num_bin)
    size_bin = np.linspace(1.0,max(big_data[i]),num_bin)
    print('num_bin: ',num_bin)

    # hist,bin = np.histogram(m250,bins= 10 * np.linspace(0.0,np.log10(max(m250)),num_bin))
    hist,bin = np.histogram(big_data[i],bins= np.linspace(1.0,max(big_data[i]),num_bin))
    hist = [(float(x) / ((0.25/((180/np.pi)**2))*((size_bin[3] - size_bin[2]) / 1e3))) for x in hist]

    for j in range(len(hist)):
    	hist[j] = np.log10(hist[j] * ((size_bin[j]/1e3)**(2.5)))

    plt.hist(bin[:-1],bin,weights=hist,histtype='step',label='%s' %(band[i]))

plt.gca().set_xscale("log")
# plt.gca().set_yscale("log")
plt.xlim((1.0,100))
plt.ylim((3.0,5.0))
plt.xlabel('S [mJy]')
plt.ylabel('log(dN/dS S$^2.5$) [Jy$^1.5$/sr$^-1$]')
plt.title('Conley Sims')
plt.legend()
plt.show()
