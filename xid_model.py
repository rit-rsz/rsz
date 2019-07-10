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

class Xid_Model():
    def __init__(self, json_dir):
        self.data = []
        for file in os.listdir(json_dir):
            if file.startswith('xid') and file.endswith('.json'):
                with open(file) as json_file:
                    datastore = json.load(json_file)
                    self.data.append(datastore)

    def plot_pixel_x_y(self, maps):
        fig, axs = plt.subplots(1,3)
        for i in range(len(self.data)):
            hdul = fits.open(maps[i]['file'])
            print(hdul[1].header['naxis1'])
            print(hdul[1].header['naxis2'])
            w = WCS(hdul[1].header)
            ra = np.array(self.data[i]['sra']) * u.deg
            dec = np.array(self.data[i]['sdec']) * u.deg
            c = SkyCoord(ra, dec)
            px, py = skycoord_to_pixel(c, w)
            axs[i].scatter(px, py)
            axs[i].set_title('posterior_model_%s_%s' % (maps[i]['name'], maps[i]['band']))
        fig.savefig('xid_image_models.png')
        plt.show()
        # print(self.data[1]['x'])
        # print(len(self.data[1]['x']))


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
    model.plot_pixel_x_y(maps)
