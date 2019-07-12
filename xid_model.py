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
from xid_test import *
import json
import matplotlib.pyplot as plt
from astropy.io import fits
from get_data import get_data
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.modeling.functional_models import Gaussian2D
from scipy.signal import convolve2d

class Xid_Model():
    def __init__(self, json_dir):
        self.data = []
        for file in os.listdir(json_dir):
            if file.startswith('xid') and file.endswith('.json'):
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
        ra = np.array(data['ra']) * u.deg
        dec = np.array(data['dec']) * u.deg
        c = SkyCoord(ra, dec)
        px, py = skycoord_to_pixel(c, w)
        plot = plt.scatter(px, py)
        plt.savefig('initial_catalog_for_xid')
        # plt.show()

    def create_psfs(self):
        self.psfs = []
        bands = [18, 25, 36]
        for i in range(len(self.data)):
            self.psfs.append([])
            for j in range(len(self.data[i]['sflux'])):
                kern = Gaussian2DKernel(bands[i] / 2.355)
                # psf = Gaussian2D(amplitude=self.data[i]['sflux'][j], x_stddev=bands[i])
                self.psfs[i].append(kern)


    def mapping_psfs(self, maps):
        gal_clusts = []
        for i in range(len(self.psfs)):
            gal_clusts.append([])
            naxis = maps[i]['astr']['NAXIS']
            gal_clust = np.empty(naxis)
            hdul = fits.open(maps[i]['file'])
            w = WCS(hdul[1].header)
            print(naxis)
            ra = np.array(self.data[i]['sra']) * u.deg
            dec = np.array(self.data[i]['sdec']) * u.deg
            c = SkyCoord(ra, dec)
            px, py = skycoord_to_pixel(c, w)
            for j in range(len(self.psfs[i])):
                psf_shape = np.asarray(self.psfs[i][j]).shape
                # print(psf_shape)
                # print(px[j]-int(psf_shape[0]/2), px[j]+int(psf_shape[0]/2), py[j]-int(psf_shape[1]/2), py[j]-int(psf_shape[1]/2))
                gal_clust = convolve2d(gal_clust[int(px[j])-int(psf_shape[0]/2):int(px[j])+int(psf_shape[0]/2), int(py[j])-int(psf_shape[1]/2):int(py[j])-int(psf_shape[1]/2)], self.psfs[i][j])
            gal_clusts.append(gal_clust)
            print(gal_clust)
            plt.show()

# def create_xid_map(priors, posterior):
#
#     """
#
#     :param priors: list of xidplus.prior classes
#     :param posterior: xidplus.posterior class
#     :return: the default xidplus Bayesian P value map plot
#     """
#     sns.set_style("white")
#     cols = ['PSW', 'PMW'. 'PLW']
#     mod_map_array = postmaps.replicated_maps(priors, posterior, posterior.samples['lp__'].size)
#     Bayes_pvals = []
#
#     cmap = sns.diverging_palette(220, 20, as_cmap=True)
#
#     hdulists = list(map(lambda prior: postmaps.make_fits_image(prior, prior.sim), priors))
#     fig = plt.figure(figsize=(10 * len(priors), 10))
#     figs = []
#     for i in range(0, len(priors)):
#         figs.append(aplpy.FITSFigure(hdulists[i][1], figure=fig, subplot=(1, len(priors), i + 1)))
#         Bayes_pvals.append(postmaps.make_Bayesian_pval_maps(priors[i], mod_map_array[i]))
#
#     for i in range(0, len(priors)):
#         figs[i].show_markers(priors[i].sra, priors[i].sdec, edgecolor='black', facecolor='black',
#                              marker='o', s=20, alpha=0.5)
#         # figs[i].tick_labels.set_xformat('dd.dd')
#         # figs[i].tick_labels.set_yformat('dd.dd')
#         figs[i]._data[
#             priors[i].sy_pix - np.min(priors[i].sy_pix) - 1, priors[i].sx_pix - np.min(priors[i].sx_pix) - 1] = \
#         Bayes_pvals[i]
#         figs[i].show_colorscale(vmin=-6, vmax=6, cmap=cmap)
#         figs[i].add_colorbar()
#         figs[i].colorbar.set_location('top')
#
#     figs[0].save('xid_all_bands_image.jpg')





if __name__ == '__main__':
    maps, err = get_data('a0370')
    model = Xid_Model('/home/vaughan/rsz/')
    model.plot_in_cat('cat_file.json', maps)
    model.plot_pixel_x_y(maps)
    model.create_psfs()
    model.mapping_psfs(maps)
