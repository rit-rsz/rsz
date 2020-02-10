import pandas as pd
import matplotlib.pyplot as plt
import time
import numpy as np
from scipy.interpolate import UnivariateSpline, interp1d
import math, sys

def ret_sides():
	a = time.time()
	data = pd.read_csv('/home/butler/SIDES_Bethermin2017.csv',usecols=[0,1,2,3,12,13,14,15,16,17])

	print('N sources : ', len(data['1.386314']))
	#retreive ra / dec pairs for each source
	ra = data['1.386314']
	dec = data['0.578664']
	z = data['0.027082']

	# retrieve 3 band fluxes
	m250 = data['0.0.5']
	m350 = data['0.0.6']
	m500 = data['0.0.7']
	# I think everything is in Jy, but leaving just in case...
	# m250 = [x * 1e3 for x in m250]
	# m350 = [x * 1e3 for x in m350]
	# m500 = [x * 1e3 for x in m500]

	b = time.time()
	print('total time: ',(b-a))

	return ra, dec, z, m250, m350, m500

''' For Plotting Histrograms
band = ['250 micron','350 micron','500 micron']
big_data = [m250,m350,m500]

for i in range(len(big_data)):
	num_bin = math.ceil(max(big_data[i])-min(big_data[i])) - 1
	# size_bin = 10 * np.linspace(0.0,np.log10(max(m250)),num_bin)
	size_bin = np.linspace(1.0,max(big_data[i]),num_bin)
	print('num_bin: ',num_bin)

	# hist,bin = np.histogram(m250,bins= 10 * np.linspace(0.0,np.log10(max(m250)),num_bin))
	hist,bin = np.histogram(big_data[i],bins= np.linspace(1.0,max(big_data[i]),num_bin))
	hist = [(float(x) / ((2.0/((180/np.pi)**2))*((size_bin[3] - size_bin[2]) / 1e3))) for x in hist]
	print(np.diff(size_bin,n=1)[0:10])
	for j in range(len(hist)):
		hist[j] = np.log10(hist[j] * ((size_bin[j]/1e3)**(2.5)))

# bins = bin[:-1] + (bin[1]-bin[0])/2
# f = UnivariateSpline(bins,hist)
# plt.scatter(bins,f(bins),color='orange')
# plt.plot(bins,f(bins),color='black')

	plt.hist(bin[:-1],bin,weights=hist,histtype='step',label='%s' %(band[i]))

plt.gca().set_xscale("log")
# y_val = interp1d(size_bin,hist)
# plt.scatter(y_val(10.0))
# plt.gca().set_yscale("log")
plt.xlim((0.1,100))
plt.ylim((3.0,5.0))
plt.xlabel('S [mJy]')
plt.ylabel('log(dN/dS S$^2.5$) [Jy$^1.5$/sr$^-1$]')
plt.title('SIDES Catalog Sims')
plt.legend()
plt.show()
'''
if __name__ == "__main__":
	ret_sides()
