import numpy as np
from astropy.io import fits
import math
#import glob
#import matplotlib
#matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

real_imagez = '/Users/vlb9398/Downloads/real_PSW.fits'
# real_image = '/Users/vlb9398/Downloads/sz_PSW.fits'
real2_image = '/Users/vlb9398/Downloads/sim_PSW.fits'
noise_image = '/Users/vlb9398/Downloads/noise_PSW.fits'
both_image = '/Users/vlb9398/Downloads/both_PSW.fits'

real_imageza = '/Users/vlb9398/Downloads/real_PMW.fits'
# real_imagea = '/Users/vlb9398/Downloads/sz_PMW.fits'
real2_imagea = '/Users/vlb9398/Downloads/sim_PMW.fits'
noise_imagea = '/Users/vlb9398/Downloads/noise_PMW.fits'
both_imagea = '/Users/vlb9398/Downloads/both_PMW.fits'

real_imagezb = '/Users/vlb9398/Downloads/real_PLW.fits'
# real_imageb = '/Users/vlb9398/Downloads/sz_PLW.fits'
real2_imageb = '/Users/vlb9398/Downloads/sim_PLW.fits'
noise_imageb = '/Users/vlb9398/Downloads/noise_PLW.fits'
both_imageb = '/Users/vlb9398/Downloads/both_PLW.fits'

bin_size = 1.0/2000.0

''' PSW '''
# real
realz = [x for x in fits.getdata(real_imagez).flatten() if not math.isnan(x)]
real_medz = np.median(realz)
fin_realz = [x - real_medz for x in realz]
hist_real,bin_real = np.histogram(fin_realz,bins= np.arange(min(fin_realz),max(fin_realz),bin_size))

# noise+sim
real = [x for x in fits.getdata(both_image).flatten() if not math.isnan(x)]
real_med = np.median(real)
fin_real = [x - real_med for x in real]
hist1,bin1 = np.histogram(fin_real,bins= np.arange(min(fin_real),max(fin_real),bin_size))

# sim
real2 = [x for x in fits.getdata(real2_image).flatten() if not math.isnan(x)]
real2_med = np.median(real2)
fin_real2 = [(x) - real2_med for x in real2]
hist3,bin3 = np.histogram(fin_real2,bins= np.arange(min(fin_real2),max(fin_real2),bin_size))

# noise
noise = [x for x in fits.getdata(noise_image).flatten() if not math.isnan(x)]
noise_med = np.median(noise)
fin_noise = [(x) - noise_med for x in noise]
hist2,bin2 = np.histogram(fin_noise,bins= np.arange(min(fin_noise),max(fin_noise),bin_size))
''' ########################################################################################## '''

''' PMW '''
# real
realza = [x for x in fits.getdata(real_imageza).flatten() if not math.isnan(x)]
real_medza = np.median(realza)
fin_realza = [x - real_medza for x in realza]
hist_reala,bin_reala = np.histogram(fin_realza,bins= np.arange(min(fin_realza),max(fin_realza),bin_size))

# noise+sim
reala = [x for x in fits.getdata(both_imagea).flatten() if not math.isnan(x)]
real_meda = np.median(reala)
fin_reala = [x - real_meda for x in reala]
hist1a,bin1a = np.histogram(fin_reala,bins= np.arange(min(fin_reala),max(fin_reala),bin_size))

# sim
real2a = [x for x in fits.getdata(real2_imagea).flatten() if not math.isnan(x)]
real2_meda = np.median(real2a)
fin_real2a = [(x) - real2_meda for x in real2a]
hist3a,bin3a = np.histogram(fin_real2a,bins= np.arange(min(fin_real2a),max(fin_real2a),bin_size))

# noise
noisea = [x for x in fits.getdata(noise_imagea).flatten() if not math.isnan(x)]
noise_meda = np.median(noisea)
fin_noisea = [(x) - noise_meda for x in noisea]
hist2a,bin2a = np.histogram(fin_noisea,bins= np.arange(min(fin_noisea),max(fin_noisea),bin_size))
''' ########################################################################################## '''

''' PLW '''
# real
realzb = [x for x in fits.getdata(real_imagezb).flatten() if not math.isnan(x)]
real_medzb = np.median(realzb)
fin_realzb = [x - real_medzb for x in realzb]
hist_realb,bin_realb = np.histogram(fin_realzb,bins= np.arange(min(fin_realzb),max(fin_realzb),bin_size))

# noise+sim
realb = [x for x in fits.getdata(both_imageb).flatten() if not math.isnan(x)]
real_medb = np.median(realb)
fin_realb = [x - real_medb for x in realb]
hist1b,bin1b = np.histogram(fin_realb,bins= np.arange(min(fin_realb),max(fin_realb),bin_size))

# sim
real2b = [x for x in fits.getdata(real2_imageb).flatten() if not math.isnan(x)]
real2_medb = np.median(real2b)
fin_real2b = [(x) - real2_medb for x in real2b]
hist3b,bin3b = np.histogram(fin_real2b,bins= np.arange(min(fin_real2b),max(fin_real2b),bin_size))

# noise
noiseb = [x for x in fits.getdata(noise_imageb).flatten() if not math.isnan(x)]
noise_medb = np.median(noiseb)
fin_noiseb = [(x) - noise_medb for x in noiseb]
hist2b,bin2b = np.histogram(fin_noiseb,bins= np.arange(min(fin_noiseb),max(fin_noiseb),bin_size))
''' ########################################################################################## '''

# Creates four polar axes, and accesses them through the returned array
fig, axes = plt.subplots(3, 3, sharey='row',sharex='row')
axes[0, 0].hist(bin3[:-1],bin3,weights=hist3,label='Sim',histtype='step')
axes[0, 0].hist(bin_real[:-1],bin_real,weights=hist_real,label='Real',histtype='step')
axes[0, 0].set_title('Lensed Sim')
axes[0, 0].set_ylabel('PSW')
axes[0, 1].hist(bin2[:-1],bin2,weights=hist2,histtype='step')
axes[0, 1].hist(bin_real[:-1],bin_real,weights=hist_real,histtype='step')
axes[0, 1].set_title('Inst Noise Only')
axes[0, 2].hist(bin1[:-1],bin1,weights=hist1,histtype='step')
axes[0, 2].hist(bin_real[:-1],bin_real,weights=hist_real,histtype='step')
axes[0, 2].set_title('Noise + Sim')

axes[1, 0].hist(bin3a[:-1],bin3a,weights=hist3a,histtype='step')
axes[1, 0].hist(bin_reala[:-1],bin_reala,weights=hist_reala,histtype='step')
axes[1, 0].set_ylabel('PMW')
axes[1, 1].hist(bin2a[:-1],bin2a,weights=hist2a,histtype='step')
axes[1, 1].hist(bin_reala[:-1],bin_reala,weights=hist_reala,histtype='step')
axes[1, 2].hist(bin1a[:-1],bin1a,weights=hist1a,histtype='step')
axes[1, 2].hist(bin_reala[:-1],bin_reala,weights=hist_reala,histtype='step')

axes[2, 0].hist(bin3b[:-1],bin3b,weights=hist3b,histtype='step')
axes[2, 0].hist(bin_realb[:-1],bin_realb,weights=hist_realb,histtype='step')
axes[2, 0].set_ylabel('PLW')
axes[2, 1].hist(bin2b[:-1],bin2b,weights=hist2b,histtype='step')
axes[2, 1].hist(bin_realb[:-1],bin_realb,weights=hist_realb,histtype='step')
axes[2, 2].hist(bin1b[:-1],bin1b,weights=hist1b,histtype='step')
axes[2, 2].hist(bin_realb[:-1],bin_realb,weights=hist_realb,histtype='step')

fig.legend(loc='upper left')
fig.text(0.5, 0.03, 'I [Jy]', ha='center')
fig.text(0.03, 0.5, 'Number/Pixel', va='center', rotation='vertical')
plt.show()
