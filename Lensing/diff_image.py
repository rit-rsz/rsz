import numpy as np
import sys, time, subprocess, os
sys.path.append('../utilities')
import config
from astropy.io import fits
import matplotlib.pyplot as plt

med_psw = config.HOME + 'Lensing/lense_template_PSW.fits'
med_pmw = config.HOME + 'Lensing/lense_template_PMW.fits'
med_plw = config.HOME + 'Lensing/lense_template_PLW.fits'

mean_psw = config.HOME + 'Lensing/lense_template_mean_PSW.fits'
mean_pmw = config.HOME + 'Lensing/lense_template_mean_PMW.fits'
mean_plw = config.HOME + 'Lensing/lense_template_mean_PLW.fits'

diff_psw = fits.getdata(med_psw) - fits.getdata(mean_psw)
diff_pmw = fits.getdata(med_pmw) - fits.getdata(mean_pmw)
diff_plw = fits.getdata(med_plw) - fits.getdata(mean_plw)

plt.imshow(diff_psw,origin = 0)
plt.colorbar()
plt.title('Diff Lensing PSW')
plt.savefig('diff_lensing_PSW.png')
plt.clf()

plt.imshow(diff_pmw,origin = 0)
plt.colorbar()
plt.title('Diff Lensing PMW')
plt.savefig('diff_lensing_PMW.png')
plt.clf()

plt.imshow(diff_plw,origin = 0)
plt.colorbar()
plt.title('Diff Lensing PLW')
plt.savefig('diff_lensing_PLW.png')
plt.clf()
