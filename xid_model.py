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

np.set_printoptions(threshold=sys.maxsize)

class Xid_Model():
    def __init__(self, json_dir, clusname):
        self.data = []
        for file in os.listdir(json_dir):
            if file.startswith('xid') and file.endswith('.json') and clusname in file:
                print(file)
                with open(file) as json_file:
                    datastore = json.load(json_file)
                    self.data.append(datastore)

    def plot_pixel_x_y(self, maps):
        # fig, axs = plt.subplots(1,3)
        for i in range(len(self.data)):
            hdul = fits.open(maps[i]['file'])
            w = WCS(hdul[1].header)
            ra = np.array(self.data[i]['sra']) * u.deg
            dec = np.array(self.data[i]['sdec']) * u.deg
            self.flux = np.array(self.data[i]['sflux'])
            c = SkyCoord(ra, dec)
            px, py = skycoord_to_pixel(c, w)
            plot = plt.scatter(px, py, c=self.flux, alpha=0.5)
            colorbar = plt.colorbar()
            colorbar.set_label('Flux')
            plt.title('XID_output_catalog_for_%s_%s' % (maps[i]['name'],maps[i]['band']))
            # axs[i].scatter(px, py, c=flux, alpha=0.5)
            # axs[i].set_title('posterior_model_%s_%s' % (maps[i]['name'], maps[i]['band']))
            # plt.set_title('posterior_model_%s_%s' % (maps[i]['name'], maps[i]['band']))
            plt.savefig('posterior_model_%s_%s' % (maps[i]['name'], maps[i]['band']))
            plt.show()
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
        initial_px = np.array(data['ra']) #* u.deg
        initial_py = np.array(data['dec']) #* u.deg
        self.cat_flux = np.array(data['flux'])
        # c = SkyCoord(ra, dec)
        # px, py = skycoord_to_pixel(c, w)
        c = pixel_to_skycoord(initial_px, initial_py, w)
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
        # print('initial', initial_px, 'final', px)

        # for i in range(len(initial_px)):
        #     print('initial:', initial_px[i], 'final', final_px[i])
        # check = []
        plt.scatter(px, py)
        plt.show()
        # plt.scatter(final_px, final_py)
        # plt.show()

        # plot = plt.scatter(ra, dec)
        plt.savefig('initial_catalog_for_xid')
        plt.show()
        # apertures = CircularAperture((ra,dec), r=4.)
        # plt.imshow(hdul[1].data, origin='lower')
        # plt.show()

    def create_psfs(self, maps):
        self.psfs = []
        bands = [18, 25, 36]
        # pixsize = [6, 25/3, 12]
        self.fluxes = []
        for i in range(len(self.data)):
        # for i in range(1):
            print(i)
            sigma = bands[i] / (sqrt(8 * log(2)))
            self.psfs.append([])
            fluxes = self.data[i]['sflux']
            self.fluxes.append(fluxes)
            hdul = fits.open(maps[i]['file'])
            w = WCS(hdul[1].header)
            ra = np.array(self.data[i]['sra']) * u.deg
            dec = np.array(self.data[i]['sdec']) * u.deg
            c = SkyCoord(ra, dec)
            px, py = skycoord_to_pixel(c, w, 1)
            # x_gen = round(bands[i] * 5 / pixsize[i])
            # if x_gen % 2 != 1:
            #     x_gen +=1
            # x_gen = int(ceil(x_gen))
            # y_gen = x_gen
            y_size = hdul[1].data.shape[0]
            x_size = hdul[1].data.shape[1]
            # if x_gen % 2 != 1:
            #     print('WARNING: npixx not odd, so PSF will not be centered')
            # if y_gen % 2 != 1:
            #     print('WARNING: npixy not odd, so psf will not be centered')
            for j in range(len(self.data[i]['sflux'])):
            # for j in range(1):
                # kern = Gaussian2DKernel(sigma, x_size=x_gen*7, y_size=y_gen*7, x_mean=px[j], y_mean=py[j])
                # kern = Gaussian2D(1, px, py, sigma, sigma)
                # kern = make_gaussian(sigma, x_gen, px[j], py[])
                # kern = np.asarray(kern)
                kern = makeGaussian(x_size=x_size, y_size=y_size, fwhm =bands[i] / maps[i]['pixsize'], center=(px[j],py[j]))
                kern = np.asarray(kern)
                # kern = kern / kern.max()
                # plt.imshow(kern)
                # plt.show()
                # kern = rebin(kern, (x_gen,y_gen))
                # print(fluxes[j])
                coefficient = fluxes[j]#fluxes[j] #/ (sqrt(pi) * bands[i] / pixsize[i]) #pi * (bands[i]/pixsize[i])**2 / (16 * log(2))#(sqrt(2*pi) * sigma / pixsize[i])
                psf = kern * coefficient
                new_tf = np.sum(psf)
                A = new_tf / fluxes[j]
                psf = psf / A
                print('Original total flux', fluxes[j])
                print('New total flux', np.sum(psf))
                # if j == self.index[0]:
                #     print('ha')
                    # plt.imshow(psf)
                    # plt.show()
                # psf = psf / psf.max()
                # psf = rebin(psf,(15, 15))
                # plt.imshow(psf)
                # plt.show()
                # psf = Gaussian2D(amplitude=self.data[i]['sflux'][j], x_stddev=bands[i])
                self.psfs[i].append(psf)


    def mapping_psfs(self, maps):
        gal_clusts = []
        print('generating mask')
        for i in range(len(self.psfs)):
            print('Starting on %s for %s elements' % (maps[i]['band'], len(self.psfs[i])))
            naxis = maps[i]['astr']['NAXIS']
            print(naxis)
            print(maps[i]['signal'].shape)
            gal_clust = np.zeros( naxis)
            print(maps[i]['file'])
            hdul = fits.open(maps[i]['file'])
            w = WCS(hdul[1].header)
            ra = np.array(self.data[i]['sra']) * u.deg
            dec = np.array(self.data[i]['sdec']) * u.deg
            c = SkyCoord(ra, dec)
            px, py = skycoord_to_pixel(c, w, 1)
            # plt.scatter(px, py)
            # c = pixel_to_skycoord(112, 162, w)
            # plt.show()
            apertures = CircularAperture((px,py), r=.4)
            plt.imshow(maps[i]['signal'].data, origin='lower')
            apertures.plot(color='red', lw=1.5, alpha=0.5)
            plt.show()
            # print(px[self.index], py[self.index])
            hdu = fits.PrimaryHDU(maps[i]['signal'].data)
            hdu = fits.HDUList([hdu])
            # hdu.writeto('original_%s_%s.fits' % (maps[i]['name'], maps[i]['band']))
            for j in range(len(self.psfs[i])):
                # print('generating mask')
                psf_shape = np.asarray(self.psfs[i][j]).shape
                # px[j] = px[j] - 1
                # py[j] = py[j] - 1
                # if int(px[j])-int(psf_shape[0]/2) >= 0 and int(px[j])+int(psf_shape[0]/2)+1 <= naxis[0] and int(py[j])-int(psf_shape[1]/2) >= 0 and int(py[j])+int(psf_shape[1]/2)+1 <= naxis[1]:
                #     x1 = floor(px[j]) - int(psf_shape[0]/2)
                #     x2 = floor(px[j]) + int(psf_shape[0]/2) + 1
                #     y1 = floor(py[j]) - int(psf_shape[1]/2)
                #     y2 = floor(py[j]) + int(psf_shape[1]/2) + 1
                # if x1 < 0:
                #     x1 = 0
                # if x2 > naxis[0]:
                #     x2 = naxis[0]
                # if y1 < 0:
                #     y1 = 0
                # if y2 > naxis[1]:
                #     y2 = naxis[1]

                gal_clust = gal_clust + self.psfs[i][j]
            # for j in range(gal_clust.shape[0]):
            #     for k in range(gal_clust.shape[1]):
            #         if gal_clust[j,k] > .46:
            #             gal_clust[j,k] = 0
            hdu = fits.PrimaryHDU(gal_clust, hdul[1].header)
            hdul = fits.HDUList([hdu])
            hdul.writeto('xid_mask_%s_%s.fits' % (maps[i]['name'], maps[i]['band']))
            print('finished generating mask for %s' % (maps[i]['band']))
            gal_clusts.append(gal_clust)
            plt.imshow(gal_clust)
            # plt.show()
                # print(j)
            # print(gal_clust)
            # plt.imshow(gal_clust)
            # plt.show()
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
        for i in range(ncols):
            if verbose:
                print('Subtracting for %s' %(maps[i]['band']))
            datafull = np.empty(maps[i]['astr']['NAXIS'])
            datasub = np.empty(maps[i]['astr']['NAXIS'])
            counter = 0
            whgd = []
            # whpl = np.where(np.isfinite(maps[i]['signal'] == False))
            for j in range(maps[i]['signal'].data.shape[0]):
                for k in range(maps[i]['signal'].data.shape[1]):
                    if np.isfinite(maps[i]['signal'].data[j,k]) == False:
                        pass
                    else:
                        whgd.append([j,k])
            # for value in whgd:
            print(models[i].shape)
            datasub = maps[i]['signal'].data - models[i]
            datafull = maps[i]['signal'].data

            whpl = []
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
         #
            # plt.imshow(datasub)
            # plt.show()

            hdu = fits.PrimaryHDU(datasub)
            hdul = fits.HDUList([hdu])
            hdul.writeto('xid_subtracted_%s_%s.fits' % (maps[i]['name'], maps[i]['band']))
         # # CONTOUR,datafull,/NODATA,$
         #         TITLE='clus_subtract_cat: Signal map for ' + $
         #         (*maps[icol]).band
         # tvimage,bytscl(datafull,min=-0.01,max=0.01,/NAN),/OVERPLOT
         #
         # CONTOUR,datasub,/NODATA,$
         #         TITLE='clus_subtract_cat: Catalog subtracted map for ' + $
         #         (*maps[icol]).band
         # tvimage,bytscl(datasub,min=-0.01,max=0.01,/NAN),/OVERPLOT

            maps[i]['scrrm'] = datasub
            maps[i]['xclean'] = datasub
        return maps

    def starfinder(self, data, fwhm):
        # Determine statistics for the starfinder application to utilize
        mean, median, std = sigma_clipped_stats(data, sigma=3.0)


        findstars = DAOStarFinder(fwhm=fwhm, threshold=1.*std)
        sources = findstars(data)
        for col in sources.colnames:
            sources[col].info.format = '%.8g'  # for consistent table output
        positions = sources['xcentroid'], sources['ycentroid']
        fluxes = sources['flux']
        # apertures = CircularAperture(positions, r=4.)
        # norm = ImageNormalize(stretch=SqrtStretch())
        # print('apertures',apertures)
        return positions, fluxes

    def plot_starfinder_flux(self, maps):
        bands = [18, 25, 36]
        pixsize = [6, 25/3, 12]
        fwhm = np.divide(bands, pixsize)
        for i in range(len(maps)):
            hdul = fits.open(maps[i]['file'])
            data = hdul[1].data
            p, f = self.starfinder(data, fwhm[i])
            x = p[0]
            y = p[1]
            apertures = CircularAperture((x,y), r=1)
            plt.imshow(data, origin='lower')
            apertures.plot()
            plt.show()
            plot = plt.scatter(x, y, c=f, alpha=0.5)
            colorbar = plt.colorbar()
            colorbar.set_label('Flux')
            plt.show()

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
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]

#
    return np.exp(-1 * ((x-x0)**2 + (y-y0)**2) / (2*sigma**2))




if __name__ == '__main__':
    # rubiks_cube()
    numpy_where_test()
    maps, err = get_data('a0370')
    model = Xid_Model('/home/vaughan/rsz/', 'a0370')
    # model.finding_index(maps)
    # print('starfinder map')
    # model.plot_starfinder_flux(maps)
    model.plot_in_cat('cat_file.json', maps)
    model.plot_pixel_x_y(maps)
    model.create_psfs(maps)
    model.find_normalization_factor(maps)
    models = model.mapping_psfs(maps)
    model.subtract_cat(maps, models)
