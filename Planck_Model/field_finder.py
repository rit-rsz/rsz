################################################################################
# NAME : ciruss_sims.py
# DATE STARTED : Julu 20, 2020
# AUTHORS : Benjamin Vaughan
# PURPOSE : To create a data structure to implement into PCAT
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import sys
sys.path.append('../utilities')
from make_map import interp_back_to_ref, create_map
from clus_get_data import *
from astropy.io import fits
import numpy as np
from math_functions import *
from planck_utilities import *
from scipy.interpolate import griddata
from astropy.wcs import WCS as world
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
import matplotlib.pyplot as plt
from astropy_healpix import HEALPix
from astropy_healpix import healpy
from astropy import units as u
from astropy.coordinates import SkyCoord, Galactic
from astropy.coordinates.representation import UnitSphericalRepresentation
from scipy.interpolate import griddata
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
import config

class field_finder():
    def __init__(self, ref_head, nu, precision, step_size, nside, nsims):
        #---- Initialize Astrometry
        self.r_pix = 3600 *  np.mean([abs(ref_head['CD1_1'] + ref_head['CD2_1']), abs(ref_head['CD2_1'] + ref_head['CD2_2'])])
        self.x_size = ref_head['naxis2']
        self.y_size = ref_head['naxis1'] #x amd y in numpy are flipped compared to fits images
        self.astrometry = ref_head

        #---- Initialize Parameters
        self.nu = nu #wavelength to calculate cirrus template at
        self.precision = precision #accuracy for how close RMS of fields have to be
        self.step_size = step_size #Gaussian step size for the field search
        self.nside = nside #this is the nside for the cut out selector not the Planck maps
        self.nsims = nsims #the number of sims to generate
        #Pull out reference Planck Field
        print('Reading in reference Planck Field.')
        PSW_I_map, ra, dec, f, b, c =  create_map(ref_head, nu=self.nu)
        self.rms_cond = round(np.sqrt(np.mean(PSW_I_map)), self.precision) #RMS conidition for 250 um
        self.sfb_comd = np.sum(PSW_I_map) #Surface Brightness Condition for 250 um
        print('Reference RMS: %.2E MJy / Sr.' % self.rms_cond)
        #----- pull in the Parameter files and store for later use.
        tau_name = config.PLANCKDATA + 'COM_CompMap_Dust-GNILC-Model-Opacity_2048_R2.01.fits'
        temp_name = config.PLANCKDATA + 'COM_CompMap_Dust-GNILC-Model-Temperature_2048_R2.01.fits'
        beta_name = config.PLANCKDATA + 'COM_CompMap_Dust-GNILC-Model-Spectral-Index_2048_R2.01.fits'

        print('Reading in Planck Parameter maps.')
        filenames = [tau_name, temp_name, beta_name]
        self.full_params = []
        for f in filenames:
            hdul = fits.open(f)
            if 'Temperature' in f:
                data = hdul[1].data.field('TEMP')
            elif 'Spectral-Index' in f:
                data = hdul[1].data.field('BETA')
            elif 'Opacity' in f:
                data = hdul[1].data.field('TAU353')
            self.full_params.append(data)

        #----- run the algorthim
        print('Searching for similar fields.')
        self.field_search()

    def conditions(self, rms, id):
        rms = round(rms, self.precision)
        if rms == self.rms_cond and id not in self.id_list:
            return True
        else:
            return False



    def pull_values(self, id_coords):
        #use id_coordinates as center values of new map.
        ra_cent  = float(id_coords.ra.to_string(decimal=True))
        dec_cent =  float(id_coords.dec.to_string(decimal=True))

        ref_head = self.astrometry #make a copy of the reference header, need to change center values.
        ref_head['CRVAL1'] = ra_cent
        ref_head['CRVAL2'] = dec_cent

        #Galactic Coordinate System
        hp = HEALPix(nside=2048, order='ring', frame='galactic') #this is with nside for the Planck Maps for the Planck Maps.

        #create a grid representing the potential map
        planck_pix = hp.pixel_resolution.to(u.arcsecond).value #pixel size of the planck parameter maps

        #the 1.5 factor is to oversize the maps
        map_deg_x = self.x_size * self.r_pix / (1.5 * 3600) #mapsize in degrees.
        map_deg_y = self.y_size * self.r_pix / (1.5 * 3600)

        npixxside =  ceil(map_deg_x * 3600 / planck_pix) #number of pixels along each axis in reference map
        npixyside =  ceil(map_deg_y * 3600 / planck_pix)

        #oversized to account for interpolation later
        ra  = np.linspace(ra_cent - map_deg_x, ra_cent + map_deg_x,   npixxside) * u.deg
        dec  = np.linspace(dec_cent - map_deg_y, dec_cent + map_deg_y,  npixyside) * u.deg
        ra_grid, dec_grid = np.meshgrid(ra, dec)

        #commit coordinate grid to a skycoord class
        coords = SkyCoord(ra=ra_grid.ravel(), dec=dec_grid.ravel(), frame='icrs')

        #pull in the parameters at these specified coordiantes
        beta = hp.interpolate_bilinear_skycoord(coords, self.full_params[2])
        T    = hp.interpolate_bilinear_skycoord(coords, self.full_params[1])
        Tau  = hp.interpolate_bilinear_skycoord(coords, self.full_params[0])

        I = calc_intensity(beta, Tau, T, nu=self.nu) * 1e20 #1e20 for conversion to Mjy

        #reshape the map to a 2D array as opposed to a 1D array.
        I_map = np.reshape(I, (npixxside, npixyside))

        #interpolate to match the reference astrometry.
        interped_map, ragrd, decgrd = interp_back_to_ref(I_map, ra_grid, dec_grid, ref_head)
        #calculate the rms and surface brightness for this realization.
        rms = np.sqrt(np.mean(interped_map))

        return rms

    def field_search(self):
        #healpix class representing a reduced resolution of the Planck Maps
        hp_class = HEALPix(nside=self.nside, order='RING', frame='Galactic')

        #calculate a starting pixel healpix ID
        id = int(abs(np.random.normal(hp_class.npix / 2., self.step_size, 1)))
        #convert this ID to a longitude / latitude and then to a RA/DEC in ICRS
        lon, lat = healpy.pix2ang(self.nside, id, lonlat=True) * u.deg
        id_coords = SkyCoord(lon, lat, representation='unitspherical').transform_to('icrs')
        #initialize an empty container for accepted field maps
        self.m_centers = np.zeros((self.nsims, 2), dtype=float)
        self.id_list = []
        self.sim   = 0


        #check to see if initial map fits the given conditions
        rms = self.pull_values(id_coords)
        if self.conditions(rms, id):
            self.id_list.append(id)
            self.m_centers[0, :] = np.array([ float(id_coords.ra.to_string(decimal=True)), float(id_coords.dec.to_string(decimal=True))])
            self.sim += 1
        #iterate to find 100 similar maps
        # for i in range(10 - len(cirrus_sim_maps)): # -len(id_list) for the case that the first realization matched the conditions.
        self.while_loop(id, rms)

        print(len(self.m_centers))
        np.save(config.PLANCKDATA + 'Planck_Models/center_values.npy', self.m_centers, allow_pickle=True)


    def while_loop(self, id, rms):
        count = 0 #initialize a counter at 0 to keep track of number of iterations.
        while not self.conditions(rms, id):
            id = int(abs(np.random.normal(id, self.step_size, 1))) #calculate another ID with a gaussian stepsize
            lon, lat = healpy.pix2ang(self.nside, id, lonlat=True) * u.deg #do the same conversion to RA/DEC in ICRS as before
            id_coords = SkyCoord(lon, lat, representation='unitspherical').transform_to('icrs')
            rms = self.pull_values(id_coords) #calculate the map model then both rms and sfbrtness
            if self.conditions(rms, id): #check the conditions and append is True
                self.m_centers[self.sim, :] = np.array([ float(id_coords.ra.to_string(decimal=True)), float(id_coords.dec.to_string(decimal=True))])
                self.id_list.append(id)
                self.sim += 1
            count += 1 #increment iteration count
            if self.sim == self.nsims:
                break
            if count > 100000: #raise an error if we have gone too far down the rabbit hole
                print(self.sim)
                raise TimeoutError('Program has reached maximum number of iterations, %s, looking for similar fields, try weakening search parameters.' % count)
                break


def SPIRE_field_generator(maps):
    centers = np.load(config.PLANCKDATA + 'Planck_Models/center_values.npy', allow_pickle=True)

    PSW, ra, dec, f, b, c = create_map(maps[0]['shead'], 1200e9)
    PMW, ra, dec, f, b, c = create_map(maps[1]['shead'], 856e9)
    PLW, ra, dec, f, b, c = create_map(maps[2]['shead'], 600e9)

    phdu = fits.PrimaryHDU()
    p250 = fits.ImageHDU(PSW / maps[0]['calfac'],  maps[0]['shead'])
    p350 = fits.ImageHDU(PMW / maps[1]['calfac'],  maps[1]['shead'])
    p500 = fits.ImageHDU(PLW / maps[2]['calfac'],  maps[2]['shead'])

    hdul = fits.HDUList([phdu, p250, p350, p500])
    hdul.writeto('sim_images/' + maps[0]['name'] + '_Planck_real.fits', overwrite=True)


    min250 = np.min(PSW) / 1.5
    max250 = np.max(PSW) * 1.5
    min350 = np.min(PMW) / 1.5
    max350 = np.max(PMW) * 1.5
    min500 = np.min(PLW) / 1.5
    max500 = np.max(PLW) * 1.5

    c = 0
    for cent in centers:
        unit = 'Jy/Beam'

        #avg_T, avg_beta, avg_tau are the same for all three bands
        PSW, ra, dec, avg_T, avg_beta, avg_tau = create_map(h250, 1200e9)
        PMW, ra, dec, avg_T, avg_beta, avg_tau = create_map(h350, 856e9)
        PLW, ra, dec, avg_T, avg_beta, avg_tau = create_map(h500, 600e9)

        PSW_head = save_header(maps[0]['shead'], avg_b, avg_T, avg_tau, 'Planck250', unit, cent)
        PMW_head = save_header(maps[1]['shead'], avg_b, avg_T, avg_tau, 'Planck350', unit, cent)
        PLW_head = save_header(maps[2]['shead'], avg_b, avg_T, avg_tau, 'Planck500', unit, cent)

        phdu = fits.PrimaryHDU()
        p250 = fits.ImageHDU(PSW / maps[0]['calfac'], h250)
        p350 = fits.ImageHDU(PMW / maps[1]['calfac'], h350)
        p500 = fits.ImageHDU(PLW / maps[2]['calfac'], h500)

        hdul = fits.HDUList([phdu, p250, p350, p500])
        hdul.writeto(config.PLANCKDATA + 'Planck_Models/'+ maps[0]['name'] + '_Planck_sim_%s.fits' % c, overwrite=True)

        fig, axs = plt.subplots(1, 3)

        im1 = axs[0].imshow(PSW, origin='lower', vmin=min250, vmax=max250)
        fig.colorbar(im1, ax=axs[0])
        im2 = axs[1].imshow(PMW, origin='lower', vmin=min350, vmax=max350)
        fig.colorbar(im2, ax=axs[1])
        im3 = axs[2].imshow(PLW, origin='lower', vmin=min500, vmax=max500)
        fig.colorbar(im3, ax=axs[2])
        plt.tight_layout()
        plt.savefig(config.PLANCKDATA + 'Planck_Models/all_3_models_%s.png' % c)
        plt.close(fig)
        plt.clf()
        c += 1



if __name__ == '__main__':
    np.random.seed(0)
    maps, err = clus_get_data('rxj1347', None)
    field_finder(maps[0]['shead'], 1200e9, 1, 5000, 256, 100)
    SPIRE_field_generator(maps)
