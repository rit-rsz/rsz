# NAME : xid_model.py
# DATE STARTED : July 9, 2019
# AUTHORS : Benjamin Vaughan
# PURPOSE : a script to create a map from the XID info.
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
# from xid_test import *
import json
import sys
import matplotlib.pyplot as plt
from astropy.io import fits
sys.path.append('../source_handling')
print(sys.path)
import numpy as np
import os
from get_data import get_data
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.modeling.functional_models import Gaussian2D
from scipy.signal import convolve as convolver
from photutils.psf import IntegratedGaussianPRF
from math import *
from photutils import CircularAperture
from astropy.wcs.utils import pixel_to_skycoord
from astropy.stats import sigma_clipped_stats
from photutils import datasets
from photutils import DAOStarFinder
from scipy.integrate import dblquad
from photutils import IRAFStarFinder
from photutils import find_peaks
np.set_printoptions(threshold=sys.maxsize)

class Xid_Model():
    def __init__(self, json_dir, clusname):
        self.data = [[],[],[]]
        for file in os.listdir(json_dir):
            if file.startswith('xid') and file.endswith('.json') and clusname in file and 'PSW' in file and 'take_2' in file:
                print(file)
                with open(json_dir + file) as json_file:
                    datastore = json.load(json_file)
                    if 'PSW' in file:
                        self.data[0] = datastore
                    elif 'PMW' in file:
                        self.data[1] = datastore
                    elif 'PLW' in file:
                        self.data[2] = datastore

    def plot_pixel_x_y(self, maps):
        # fig, axs = plt.subplots(1,3)
        for i in range(1):
            hdul = fits.open(maps[1]['file'])
            w = WCS(hdul[1].header)
            ra = np.array(self.data[i]['sra']) * u.deg
            dec = np.array(self.data[i]['sdec']) * u.deg
            self.flux = np.array(self.data[i]['sflux'])
            c = SkyCoord(ra, dec)
            px, py = skycoord_to_pixel(c, w)
            plot = plt.scatter(px, py, c=self.flux, alpha=0.5)
            colorbar = plt.colorbar()
            colorbar.set_label('Flux')
            plt.title('XID_output_catalog_for_%s_%s' % (maps[0]['name'],maps[0]['band']))
            # axs[i].scatter(px, py, c=flux, alpha=0.5)
            # axs[i].set_title('posterior_model_%s_%s' % (maps[i]['name'], maps[i]['band']))
            # plt.set_title('posterior_model_%s_%s' % (maps[i]['name'], maps[i]['band']))
            # plt.savefig('posterior_model_%s_%s' % (maps[i]['name'], maps[i]['band']))
            # plt.show()
            # plt.show()
        # print(self.data[1]['x'])
        # print(len(self.data[1]['x']))

    def plot_in_cat(self, filename, maps):
        for i in range(len(maps)):
            if 'PSW' in maps[i]['band']:
                index = i
        hdul = fits.open(maps[index]['file'])
        w = WCS(hdul[1].header)
        with open(filename) as json_file:
            data = json.load(json_file)
        mean, median, std = sigma_clipped_stats(hdul[1].data, sigma=3.0)
        initial_px = np.array(data['ra']) #* u.deg
        initial_py = np.array(data['dec']) #* u.deg
        self.cat_flux = np.array(data['flux']) * std
        # c = SkyCoord(ra, dec)
        # px, py = skycoord_to_pixel(c, w)
        c = pixel_to_skycoor(initial_px, initial_py, w)
        final_px, final_py = skycoord_to_pixel(c, w)


        a = []
        d = []
        coordinates = c.to_string('decimal')
        for i in range(len(coordinates)):
            split_coords = coordinates[i].split(' ')
            a.append(float(split_coords[0]))
            d.append(float(split_coords[1]))

        c = SkyCoord(a * u.deg, d * u.deg)
        px, py = skycoord_to_pixel(c, w)
        apertures = CircularAperture((px,py), r=.4)
        plt.imshow(hdul[1].data, origin='lower')
        apertures.plot()
        plt.show()
        plt.scatter(px, py)
        plt.show()
        plt.savefig('initial_catalog_for_xid')
        plt.show()

    def create_psfs(self, maps):
        self.psfs = []
        bands = [18, 36, 25]
        # pixsize = [6, 25/3, 12]
        # self.fluxes = []
        # for i in range(3):
        for i in range(1):
            fwhm = bands[i] / maps[i]['pixsize']
            sigma = fwhm / (sqrt(8 * log(2)))
            print('fwhm in arcsec', bands[i])
            print('pixsize in arcsec / pixel', maps[i]['pixsize'])
            print('fwhm in pixels', fwhm)
            print('sigma', sigma)
            self.psfs.append([])
            fluxes = self.data[i]['sflux']
            hdul = fits.open(maps[1]['file'])
            w = WCS(hdul[1].header)
            print(self.data[i]['sra'])
            ra = np.array(self.data[i]['sra']) * u.deg
            print(ra)
            dec = np.array(self.data[i]['sdec']) * u.deg
            c = SkyCoord(ra, dec)
            px, py = skycoord_to_pixel(c, w, 1)
            print(bands[1], maps[1]['file'])
            # fluxes = np.array(fluxes)
            y_size = hdul[1].data.shape[0]
            x_size = hdul[1].data.shape[1]
            A = integrate(x_size, y_size, sigma)[0]
            # print(A[0])
            for j in range(len(fluxes)):
                kern = makeGaussian(x_size=x_size, y_size=y_size, fwhm =bands[i] / maps[i]['pixsize'], center=(px[j], py[j]))
                kern = np.asarray(kern)
                kern = kern / np.max(kern)
                coefficient = fluxes[j]
                psf = kern * coefficient
                self.psfs[i].append(psf)


    def mapping_psfs(self, maps):
        gal_clusts = []
        print('generating mask')
        for i in range(1):
            print('Starting on %s for %s elements' % (maps[1]['band'], len(self.psfs[i])))
            naxis = maps[1]['astr']['NAXIS']
            print(naxis)
            print(maps[1]['signal'].shape)
            gal_clust = np.zeros( naxis)
            print(maps[1]['file'])
            hdul = fits.open(maps[i]['file'])
            w = WCS(hdul[1].header)
            ra = np.array(self.data[i]['sra']) * u.deg
            dec = np.array(self.data[i]['sdec']) * u.deg
            c = SkyCoord(ra, dec)
            px, py = skycoord_to_pixel(c, w, 1)
            apertures = CircularAperture((px,py), r=.4)
            hdu = fits.PrimaryHDU(maps[i]['signal'].data)
            hdu = fits.HDUList([hdu])
            for j in range(len(self.psfs[i])):
                psf_shape = np.asarray(self.psfs[i][j]).shape
                gal_clust = gal_clust + self.psfs[i][j]
            gal_clusts.append(gal_clust)
            hdu = fits.PrimaryHDU(gal_clust, hdul[1].header)
            hdul = fits.HDUList([hdu])
            hdul.writeto('starfinder_total_mask_%s_%s.fits' % (maps[1]['name'], maps[1]['band']))
            print('finished generating mask for %s' % (maps[1]['band']))

        return gal_clusts

    def finding_index(self, maps):
        signal = maps[0]['signal'].data
        print(signal[165, 112])
        print(signal[170, 108])
        print(signal[164, 112])
        print(signal[182, 186])

    def subtract_cat(self, maps, models, verbose=1):
        err = False
        ncols = len(maps)
        print(len(models))
        print(len(models[0]))
        for i in range(1):
            if verbose:
                print('Subtracting for %s' %(maps[1]['band']))
            datafull = np.empty(maps[1]['astr']['NAXIS'])
            datasub = np.empty(maps[1]['astr']['NAXIS'])
            counter = 0
            whgd = []
            # whpl = np.where(np.isfinite(maps[i]['signal'] == False))
            for j in range(maps[1]['signal'].data.shape[0]):
                for k in range(maps[1]['signal'].data.shape[1]):
                    if np.isfinite(maps[1]['signal'].data[j,k]) == False:
                        pass
                    else:
                        whgd.append([j,k])
            # for value in whgd:
            print(models[i].shape)
            datasub = maps[1]['signal'].data - models[i]
            datafull = maps[1]['signal'].data
            # for j in range(maps[i]['mask'].shape[0]):    # hdu = fits.PrimaryHDU(gal_clust)
            # hdul = fits.HDUList([hdu])
            # hdul.writeto('mask_%s_%s.fits' % (maps[i]['name'], maps[i]['band']))
            # print('finished generating mask for %s' % (maps[i]['band']))
            # gal_clusts.append(gal_clust)
            #     for k in range(maps[i]['mask'].shape[1]):
            #         if maps[i]['mask'][j,k] == 1:
            #             whpl.append([j,k])
            # # whpl = np.where(maps[i]['mask'] == 1)
            # for value in whpl:
            #     datasub[value] = np.nan
         #            # self.fluxes.append(fluxes)

            # plt.imshow(datasub)
            # plt.show()

            hdu = fits.PrimaryHDU(datasub)
            hdul = fits.HDUList([hdu])
            hdul.writeto('starfinder_total_subtracted_%s_%s.fits' % (maps[1]['name'], maps[1]['band']))
         # # CONTOUR,datafull,/NODATA,$
         #         TITLE='clus_subtract_cat: Signal map for ' + $
         #         (*maps[icol]).band
         # tvimage,bytscl(datafull,min=-0.01,max=0.01,/NAN),/OVERPLOT
         #
         # CONTOUR,datasub,/NODATA,$
         #         TITLE='clus_subtract_cat: Catalog subtracted map for ' + $
         #         (*maps[icol]).b, len(x2_ind)and
         # tvimage,bytscl(datasub,min=-0.01,max=0.01,/NAN),/OVERPLOT

            maps[i]['scrrm'] = datasub
            maps[i]['xclean'] = datasub
        return maps

    def starfinder(self, data, fwhm):
        # Determine statistics for the starfinder application to utilize
        mean, median, std = sigma_clipped_stats(data, sigma=3.0)
        print(std) #0.004261477072239812
        findstars = DAOStarFinder(fwhm=fwhm, threshold=1.*std)
        sources = findstars(data)
        for col in sources.colnames:
            sources[col].info.format = '%.8g'  # for consistent table output
        positions = sources['xcentroid'], sources['ycentroid']
        fluxes = sources['flux']
        # apertures = CircularAperture(positions, r=4.)
        # norm = ImageNormalize(stretch=SqrtStretch())
        # print('apertures',apertures)
        print(std)
        print(fwhm)
        return positions, fluxes, std

    def plot_starfinder_flux(self, maps):
        self.flux = []
        bands = [25, 25, 36]
        pixsize = [6, 25/3, 12]
        fwhm = np.divide(bands, pixsize)
        for i in range(1):
            # self.flux.append([])
            hdul = fits.open('/data/mercado/SPIRE/hermes_clusters/a0370_PSW_nr_1.fits')
            print(maps[0]['file'], 'HE SCREAM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11')
            data = hdul[1].data
            p, f, std = self.starfinder(data, fwhm[i])
            x = p[0]
            y = p[1]
            f = f * 1. * std
            apertures = CircularAperture((x,y), r=1)
            plt.imshow(data, origin='lower')
            apertures.plot()
            # plt.show()
            print('THIS IS THE STARFINDER MAP BEFORE PASSSING IT TO THE NEW FUNCTION !!!!!!!!')
            plot = plt.scatter(x, y, c=f, alpha=0.5)
            colorbar = plt.colorbar()
            colorbar.set_label('Flux')
            self.sf_flux = f
            self.sf_x = x
            self.sf_y = y
            print('Length from starfinder', len(self.sf_x))
            # plt.show()

    def plot_IRAFstarfinder(self, maps):
        data = fits.open(maps[0]['file'])[1].data
        print(maps[1]['file'])
        mean, mediam, std = sigma_clipped_stats(data, sigma=3.0)
        starfinder = IRAFStarFinder(threshold=.001*std, fwhm=3.0)
        table = find_peaks(data, 2*std)
        print(table.colnames)
        for col in table.colnames:
            table[col].info.format = '%.8g'
        positions = (table['x_peak'], table['y_peak'])
        fluxes = table['peak_value']
        print(len(fluxes))
        plt.scatter(positions[0], positions[1], c=fluxes, alpha=0.5)
        colorbar = plt.colorbar()
        # plt.show()
        print('length of find_peaks', len(fluxes))
        self.peak_fluxes = fluxes
        self.peak_p = positions

    def find_normalization_factor(self, maps):
        for i in range(len(maps)):
            hdul = fits.open(maps[i]['file'])
            w = WCS(hdul[1].header)
            ra = np.array(self.data[i]['sra']) * u.deg
            dec = np.array(self.data[i]['sdec']) * u.deg
            c = SkyCoord(ra, dec)
            px, py = skycoord_to_pixel(c, w, 1)
            new_f = 0
            for j in range(len(self.fluxes[i])):
                new_f += self.fluxes[i][j] / hdul[1].data[int(py[j])-1, int(px[j])-1]
            avg = new_f / len(self.fluxes[i])
            print(avg)


    def bubble_sort(self, arr):
        n = len(arr[0])
        for i in range(n):
            for j in range(0, n-i-1):
                if arr[0][j] > arr[0][j+1]:
                    holder = arr[0][j]
                    arr[0][j] = arr[0][j+1]
                    arr[0][j+1] = arr[0][j]
                    holder = arr[1][j]
                    arr[1][j] = arr[1][j+1]
                    arr[1][j+1] = arr[1][j]
        return np.asarray(arr)


    def create_PSW_csv(self):

        # \n

        #HeDam 1 data
        f_obj = open('cesam_all_scat250_dr2_catalog_1564065377.csv')

        HeDam1_RA = []
        HeDam1_Dec = []
        HeDam1_Flux = []
        for line in f_obj:
            line = line.split(',')
            if line[2] == 'RA':
                pass
            else:
                HeDam1_RA.append(float(line[2]))
                HeDam1_Dec.append(float(line[3]))
                HeDam1_Flux.append(float(line[4]))

        HeDam1_RA = np.asarray(HeDam1_RA)
        HeDam1_Dec = np.asarray(HeDam1_Dec)
        HeDam1_Flux = np.asarray(HeDam1_Flux)

        # HeDam1_data = np.load('a2218_250.npy')
        # HeDam1_RA = HeDam1_data[0]
        # HeDam1_Dec = HeDam1_data[1]
        # HeDam1_Flux = HeDam1_data[2]
        # HeDam1 = fits.open('/home/vaughan/rsz/fits_files/CS-Abell-370_SCAT250_DR2.fits')
        # HeDam1_head = HeDam1[1].header
        # HeDam1_data = HeDam1[1].data
        # HeDam1_RA = HeDam1_data.field('RA')
        # HeDam1_Dec = HeDam1_data.field('Dec')
        # HeDam1_Flux = HeDam1_data.field('F250')
        # plt.scatter(HeDam1_RA, HeDam1_Dec, c=HeDam1_Flux, alpha=0.5)
        # plt.show()

        print(HeDam1_RA.shape)
        #HeDam2 data
        # HeDam2 = fits.open('/home/vaughan/rsz/fits_files/CS-Abell-370_SCAT250SXT_DR2.fits')
        # HeDam2_head = HeDam2[1].header
        # HeDam2_data = HeDam2[1].data
        # HeDam2_RA = HeDam2_data.field('RA')
        # HeDam2_Dec = HeDam2_data.field('Dec')
        # HeDam2_Flux = HeDam2_data.field('Flux')

        #peaks from find_peaks
        peak_x = self.peak_p[0]
        peak_y = self.peak_p[1]
        peak_f = self.peak_fluxes

        plt.scatter(peak_x, peak_y, c=peak_f, alpha=0.5)
        # plt.show()

        #fluxes from starfinder
        sf_x = self.sf_x
        sf_y = self.sf_y
        sf_f = self.sf_flux

        print('THIS IS THE STARFINDER MAP AFTER PASSING IT TO THE NEW FUNCTION !!!!')
        plt.scatter(sf_x, sf_y, c=sf_f, alpha=0.5)
        # plt.show()

        #XID
        XID_data = self.data[0]
        XID_RA = XID_data['sra'] * u.deg
        XID_Dec = XID_data['sdec'] * u.deg
        XID_Flux = XID_data['sflux']
        XID_head = fits.open('/data/mercado/SPIRE/hermes_clusters/a0370_PSW_nr_1.fits')[1].header
        w = WCS(XID_head)
        c = SkyCoord(XID_RA, XID_Dec)
        px, py = skycoord_to_pixel(c, w, 1)

        print('num of elements in xid', len(px))

        plt.scatter(px, py, c=XID_Flux, alpha=0.5)
        # plt.show()

        print(len(XID_RA))
        #Time to interpolate!!!!!!!!!

        XID_RA = XID_RA /  u.deg
        XID_Dec = XID_Dec / u.deg

        XID_RA = np.asarray(XID_RA)
        XID_Dec = np.asarray(XID_Dec)

        min_ra = .0064
        min_dec = .0064
        min_px1 = 2.5
        min_py1 = 2.5
        min_px = 2.32
        min_py = 2.32

        print('interpolating')

        holder = 100
        sf_ind = []
        for i in range(len(XID_Flux)):
            holder = 100
            for j in range(len(sf_x)):
                d = sqrt((px[i]-sf_x[j])**2 + (py[i]-sf_y[j])**2)
                if d < holder:
                    holder = d
                    index = np.where(sf_x == sf_x[j])[0][0]
            sf_ind.append(index)

        print('length of sf_ind', len(sf_ind))

        holder = 100
        peak_ind = []
        for i in range(len(XID_Flux)):
            holder = 100
            for j in range(len(peak_x)):
                d = sqrt((px[i]-peak_x[j])**2 + (py[i]-peak_y[j])**2)
                if d < holder:
                    holder = d
                    index = np.where(peak_x == peak_x[j])[0][0]
            peak_ind.append(index)

        print('num elements after interpolation', len(sf_ind))

        holder = 100
        h1_ind = []
        for i in range(len(XID_Flux)):
            holder = 100
            for j in range(len(HeDam1_RA)):
                d = sqrt((XID_RA[i]-HeDam1_RA[j])**2 + (XID_Dec[i]-HeDam1_Dec[j])**2)
                if d < holder:
                    holder = d
                    index = np.where(HeDam1_RA == HeDam1_RA[j])[0][0]
            h1_ind.append(index)

        print('num elements after interpolation', len(peak_ind))


        print('sorting arrays')

        h1_ind = np.asarray(h1_ind)
        HeDam1_Flux = HeDam1_Flux[h1_ind] / 1000
        sf_f = sf_f[sf_ind]
        peak_f = peak_f[peak_ind]

        ind = np.nonzero(HeDam1_Flux)[0]
        print(ind)

        n_HeDam1_f = HeDam1_Flux[ind]
        n_XID_f = np.asarray(XID_Flux)
        n_XID_f = n_XID_f[ind]

        percent_err = abs(n_XID_f - n_HeDam1_f) / n_HeDam1_f
        mean = np.mean(percent_err) * 100
        """
        print('h1')
        x1_ind, H1_ind = self.index_finder(XID_RA, XID_Dec, HeDam1_RA, HeDam1_Dec, min_ra, min_dec)
        print('sf')
        x3_ind, sf_ind = self.index_finder(px, py, sf_x, sf_y, min_px1, min_py1)
        print('peak')
        x4_ind, peak_ind = self.index_finder(px, py, peak_x, peak_y, min_px, min_py)
        print(len(x1_ind), len(x3_ind), len(x4_ind))

        print(x1_ind)
        print(x3_ind)
        print(x4_ind)        print(':)', len(HeDam1_RA))


        HeDam1_Flux = HeDam1_Flux[H1_ind] / 1000
        sf_f = sf_f[sf_ind]
        peak_f = peak_f[peak_ind]
        XID_Flux = np.asarray(XID_Flux)[x1_ind]

        """

        f_obj = open('flux_catalog.txt', 'w')

        f_obj.write('Catalog of Fluxes for Abell 2218 using /data/mercado/SPIRE/hermes_clusters/a2218_PSW_nr_1.fits \n' )
        f_obj.write('Mean percent error of XID flux with HeDam Flux is: %s \n' % (mean))
        f_obj.write('Source Num     Pixel x          Pixel y     HeDam Flux         XID flux     DAO flux     Peak Flux \n ')

        for i in range(len(XID_Flux)):
            str =   '%s              %.5f       %.5f       %.5f          %.5f        %.5f        %.5f  \n' % (i, px[i], py[i], HeDam1_Flux[i], XID_Flux[i],  sf_f[i], peak_f[i])
            f_obj.write(str)
        f_obj.close()



    def index_finder(self, XID_RA, XID_Dec, other_RA, other_dec, x_off, y_off):
        other_index = []
        x_index = []
        for RA, Dec in zip(XID_RA, XID_Dec):
            for r,d in zip(other_RA, other_dec):
                if RA + x_off >= r and RA - x_off <= r and Dec + y_off >= d and Dec - y_off <= d: # and Dec + y_off >= d and Dec + y_off <= d:
                    index = np.where(other_RA == r)[0][0]
                    other_index.append(index)
                    index = np.where(XID_RA == RA)[0][0]
                    x_index.append(index)
        return x_index, other_index


def rebin(a, new_shape):
    shape = a.shape
    M = int(shape[0])
    N = int(shape[1])
    m, n = new_shape
    if m<M:
        return a.reshape((m,int(M/m),n,int(N/n))).mean(3).mean(1)
    else:
        return np.repeat(np.repeat(a, m/M, axis=0), n/N, axis=1)

def rubiks_cube():
    rubiks_arr = np.array([[1 ,2 ,3 ,4 ,5 ],
                           [6 ,7 ,8 ,9 ,10],
                           [11,12,13,14,15],
                           [16,17,18,19,20],
                           [21,22,23,24,25]])
    print('index of 3,4', rubiks_arr[3,4])
    plt.imshow(rubiks_arr)
    plt.show()
    hdu = fits.PrimaryHDU(rubiks_arr)
    hdul = fits.HDUList([hdu])
    hdul.writeto('rubiks.fits')

def numpy_where_test():
    arr = np.array([[2,3,67,3,4],
                    [5,6,6,3,6]])
    index = np.where(arr == 67)
    print(index)
    print(index[0], index[1])
#
# def rotate_origin_only(xy, radians):
#     """Only rotate a point around the origin (0, 0)."""
#     x, y = xy
#     xx = x * math.cos(radians) + y * math.sin(radians)
#     yy = -x * math.sin(radians) + y * math.cos(radians)
#
# return xx, yy

def make_gaussian(sigma, size, x_mean, y_mean):
    pass

def makeGaussian(x_size, y_size, fwhm = 3, center=None):
    """ Make a square gaussian kernel.

    size is the length of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """

    sigma = fwhm / 2.355
    x = np.arange(0, x_size, 1, float)
    y = np.arange(0, y_size, 1, float)
    y = y[:,np.newaxis]

    if center is None:
        x0 = x_size // 2
        y0 = y_size // 2
    else:
        x0 = center[0]
        y0 = center[1]

    sigma = fwhm / 2.355
#
    return np.exp(-1 * ((x-x0)**2 + (y-y0)**2) / (2*sigma**2))

def integrate(x_size, y_size, sigma):
    x0 = x_size // 2
    y0 = y_size // 2
    area = dblquad(lambda x, y: np.exp(-1 * ((x-x0)**2 + (y-y0)**2) / (2*sigma**2)), 0, x_size, lambda x: 0, y_size)
    return area

def noise_map():
    im = fits.open('/home/vaughan/rsz/fits_files/original_a0370_PSW.fits')
    data = im[0].data
    map = np.empty(data.shape)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if data[i,j] >= 0:
                data[i,j] = data[i,j] * 1000
                map[i,j] = sqrt(data[i,j])
            elif data[i,j] < 0:
                data[i,j] = data[i,j] * 1000
                map[i,j] = -1 * sqrt(-1 * data[i,j])
            else:
                map[i,j] = np.nan

    hdu = fits.PrimaryHDU(map)
    hdul = fits.HDUList([hdu])
    hdul.writeto('noise_test.fits')


def variance_or_noise():
    fits_file = fits.open('/data/mercado/SPIRE/hermes_clusters/a0370_PSW_nr_1.fits')
    data = fits_file[1].data
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)

def subtract_models(model1, model2):
    residual = model1 - model2
    hdu = fits.PrimaryHDU(residual)
    hdul = fits.HDUList([hdu])
    hdul.writeto('idl-findpeaks.fits')

if __name__ == '__main__':
    # rubiks_cube()

    # variance_or_noise()
    # numpy_where_test()
    # g = makeGaussian(105, 105, fwhm=18)
    # plt.imshow(g)
    # plt.show()
    maps, err = get_data('a0370')
    print(maps[0]['file'])
    # # # noise_map()
    model = Xid_Model('/home/vaughan/rsz/json_files/', 'a0370')
    model.plot_IRAFstarfinder(maps)
    # # # # model.finding_index(maps)
    # # # # print('starfinder map')
    model.plot_starfinder_flux(maps)
    # # # # model.plot_in_cat('cat_file.json', maps)
    # # model.plot_pixel_x_y(maps)
    # model.create_psfs(maps)
    # # # # model.find_normalization_factor(maps)
    # models = model.mapping_psfs(maps)
    # model.subtract_cat(maps, models)
    model.create_PSW_csv()

    # fits1 = fits.open('/home/vaughan/rsz/fits_files/idl_subtracted_a0370_PSW.fits')
    # fits2 = fits.open('starfinder_total_subtracted_a0370_PSW.fits')
    # model1 = fits1[0].data
    # # for i in range(model1.shape[0]):
    # #     for j in range(model1.shape[1]):
    # #         if model1[i,j] >= 0:
    # #             model1[i,j] = model1[i,j] / 1000
    # #             # map[i,j] = sqrt(data[i,j])
    # #         elif model1[i,j] < 0:
    # #             model1[i,j] = model1[i,j] / 1000
    #             # map[i,j] = -1 * sqrt(-1 * data[i,j])
    # model2 = fits2[0].data
    # subtract_models(model1, model2)
    # # models = [[],[],[]]
    # for file in os.listdir('/home/vaughan/rsz'):
    #     if 'xid_model' in file and '.fits' in file:
    #         hdul = fits.open(file)
    #         data = hdul[0].data
    #         if 'PSW' in file:
    #             models[0] = data
    #         elif 'PMW' in file:
    #             models[1] = data
    #         elif 'PLW' in file:
    #             models[2] = data
    # model.subtract_cat(maps, models)
