from __future__ import print_function, absolute_import

import numpy as np
import math, sys
import astropy.io.fits as fits
import matplotlib.pyplot as plt
sys.path.append('/home/butler/rsz/new_bethermin/')

try:
    from astropy.convolution import convolve
    from astropy.convolution import Gaussian2DKernel
except ImportError:
    raise Exception("You have an old (pre 0.3) version of astropy")


""" Generates simulated maps"""

# __all__ = ["genmap_gauss", "get_gauss_beam"]

import gencat


def get_gauss_beam(fwhm, pixscale, nfwhm=5.0, oversamp=1):
    """ Generate Gaussian kernel

    Parameters
    ----------
    fwhm: float
      FWHM of the Gaussian beam.

    pixscale: float
      Pixel scale, in same units as FWHM.

    nfwhm: float
      Number of fwhm (approximately) of each dimension of the output beam.

    oversamp: int
      Odd integer giving the oversampling to use when constructing the
      beam.  The beam is generated in pixscale / oversamp size pixels,
      then rebinned to pixscale.

    Notes
    -----
      The beam is normalized by having a value of 1 in the center.
      If oversampling is used, the returned array will be the sum over
      the neighborhood of this maximum, so will not be one.
    """

    if fwhm <= 0:
        raise ValueError("Invalid (negative) FWHM")
    if pixscale <= 0:
        raise ValueError("Invalid (negative) pixel scale")
    if nfwhm <= 0.0:
        raise ValueError("Invalid (non-positive) nfwhm")
    if fwhm / pixscale < 2.5:
        raise ValueError("Insufficiently well sampled beam")
    if oversamp < 1:
        raise ValueError("Invalid (<1) oversampling")

    retext = round(fwhm * nfwhm / pixscale)
    if retext % 2 == 0:
        retext += 1

    bmsigma = fwhm / math.sqrt(8 * math.log(2))

    beam = Gaussian2DKernel(bmsigma / pixscale, x_size=retext,
                            y_size=retext, mode='oversample',
                            factor=oversamp)
    beam *= 1.0 / beam.array.max()
    return beam


class genmap_gauss :
    """ Generates simulated maps from the Bethermin et al. 2012 model
    using a Gaussian beam"""

    def __init__(self, wave=[250.0, 350, 500], pixsize=[6.0, 8.33333, 12.0],
                 fwhm=[17.6, 23.9, 35.2], nfwhm=5.0, bmoversamp=5,
                 gensize=150000, truthtable=True, log10Mb=11.2,
                 alpha=1.3, log10Mmin=8.5, log10Mmax=12.75, ninterpm=2000,
                 zmin=0.1, zmax=10.0, Om0=0.315, H0=67.7, phib0=-3.02,
                 gamma_sfmf=0.4, ninterpz=1000, rsb0=0.012, gammasb=1.0,
                 zsb=1.0, logsSFRM0=-10.2, betaMS=-0.2, zevo=2.5,
                 gammams=3.0, bsb=0.6, sigmams=0.15, sigmasb=0.2,
                 mnU_MS0=4.0, gammaU_MS0=1.3, z_UMS=2.0, mnU_SB0=35.0,
                 gammaU_SB0=0.4, z_USB=3.1, scatU=0.2, ninterpdl=200):
        """ Initializer.

        Parameters
        ----------
        wave: ndarray
          Wavelengths to generate maps at, in microns.

        pixsize: ndarray
          Pixel sizes of output maps, in arcsec.  The code may
          perform somewhat better if the finest pixel scale is first.

        fwhm: ndarray
          Beam FWHM values, in arcsec.

        nfwhm: float
          How far out, in units of FWHM, to generate the beams

        gensize: int
          In order to try to save memory, sources are added to the maps
          in chunks of this size.  If set to 0, all the sources are
          generated at once.

        truthtable: bool
          If set to true, then the truth table is also returned (giving
          the source fluxes and positions).  If this is set, then
          gensize is set to 0.

        log10Mb: float
          Log10 of Mb, in solar masses

        alpha: float
          Power-law slope of low mass distribution

        log10Mmin: float
          Log10 minimum mass to generate, in solar masses

        log10Mmax: float
          Log10 maximum mass to generate, in solar masses

        ninterpm: float
          Number of mass interpolation samples to use

        zmin: float
          Minimum z generated

        zmax: float
          Maximum z generated

        Om0: float
          Density parameter of matter

        H0: float
          Hubble constant in km / s / Mpc

        phib0: float
          log 10 number density at SFMF break in comoving Mpc^-3

        gamma_sfmf: float
          Evolution of the number density of the SFMF at z > 1

        ninterpz: int
          Number of z interpolation samples to use

        rsb0: float
          Relative amplitude of SB distribution to MS one.

        gammasb: float
          Redshift evolution in SB relative amplitude.

        zsb: float
          Redshift where SB relative amplitude stops evolving

        logsSFRM0: float
          Base Log 10 specific star formation rate.

        betaMs: float
          Slope of sFRM-M relation

        zevo: float
          Redshift where MS normalization stops evolving

        gammams: float
          Redshift evolution power law for sSFR

        bsb: float
          Boost in sSFR for starbursts, in dex

        sigmams: float
          Width of MS log-normal distribution

        sigmasb: float
          Width of SB log-normal distribution

        mnU_MS0: float
          Average ultraviolet intensity field in MS galaxies at z=0

        gammaU_MS0: float
          Evolution in <U> for MS

        z_UMS: float
          Redshift where <U> stops evolving for MS galaxies

        mnU_SB0: float
          Average ultraviolet intensity field in SB galaxies at z=0

        gammaU_SB0: float
          Evolution in <U> for SB

        z_USB: float
          Redshift where <U> stops evolving for SB galaxies

        scatU: float
          Scatter in <U>, in dex

        ninterpdl: float
          Number of interpolation points in the luminosity distance.
        """

        from numbers import Number

        if isinstance(wave, Number):
            self._wave = np.asarray([wave], dtype=np.float32)
        else:
            self._wave = np.asarray(wave, dtype=np.float32)
        if self._wave.min() <= 0:
            raise ValueError("Non-positive wavelengths not supported")
        self._nbands = len(self._wave)

        if isinstance(pixsize, Number):
            self._pixsize = pixsize * np.ones_like(self._wave,
                                                   dtype=np.float32)
        else:
            if len(pixsize) != self._nbands:
                if len(pixsize) == 1:
                    self._pixsize = np.asarray(pixsize[0],
                                               dtype=np.float32) * \
                        np.ones_like(self._wave)
                else:
                    raise ValueError("Number of pixel sizes doesn't "
                                     "match number of wavelengths")
            else:
                self._pixsize = np.asarray(pixsize, dtype=np.float32)
        if self._pixsize.min() <= 0:
            raise ValueError("Invalid (negative) pixel size")

        if isinstance(fwhm, Number):
            self._fwhm = fwhm * np.ones_like(self._wave)
        else:
            if len(fwhm) != self._nbands:
                if len(fwhm) == 1:
                    self._fwhm = np.asarray(fwhm[0], dtype=np.float32) *\
                        np.ones_like(self._wave)
                else:
                    raise ValueError("Number of FWHM doesn't match number "
                                     "of wavelengths")
            else:
                self._fwhm = np.asarray(fwhm, dtype=np.float32)
        if self._fwhm.min() <= 0:
            raise ValueError("Invalid (negative) FWHM")
        if (self._fwhm / self._pixsize).min() < 1:
            raise ValueError("Some FWHM not properly sampled by pixel size")

        self._nfwhm = float(nfwhm)
        if self._nfwhm <= 0:
            raise ValueError("Invalid (non-positive) nfwhm")

        # Control how sources are generated
        if truthtable:
            self._gensize = 0
            self._returntruth = True
        else:
            # If gensize is zero, we do the full set at once
            self._gensize = int(gensize)
            self._returntruth = False
            if self._gensize < 0:
                raise ValueError("Invalid (negative) gensize")

        self._bmoversamp = int(bmoversamp)
        if self._bmoversamp < 1:
            raise ValueError("Invalid (<1) beam oversampling")
        if self._bmoversamp % 2 == 0:
            errstr = "Invalid (even) beam oversampling {:d}"
            raise ValueError(errstr.format(self._bmoversamp))

        # Set up catalog generator
        self._gencat_stuff = gencat.Gencat(log10Mb, alpha, log10Mmin, log10Mmax,
                              ninterpm, zmin, zmax, Om0, H0, phib0,
                              gamma_sfmf, ninterpz, rsb0, gammasb, zsb,
                              logsSFRM0, betaMS, zevo, gammams, bsb,
                              sigmams, sigmasb, mnU_MS0, gammaU_MS0,
                              z_UMS, mnU_SB0, gammaU_SB0, z_USB,
                              scatU, ninterpdl)
        self._npersr = self._gencat_stuff.npersr

    @property
    def npersr(self):
        return self._npersr

    @property
    def nbands(self):
        return self._nbands

    @property
    def wave(self):
        return self._wave

    @property
    def pixsize(self):
        return self._pixsize

    @property
    def fwhm(self):
        return self._fwhm

    @property
    def bmoversamp(self):
        return self._bmoversamp

    @property
    def gensize(self):
        return self._gensize

    def get_gauss_beam(self, idx):
        """ Gets the beam for the specified index"""
        return get_gauss_beam(self._fwhm[idx], self._pixsize[idx],
                              self._nfwhm, self._bmoversamp)

    def generate(self, area, sigma=None, verbose=False):
        """ Generates simulated maps.

        Parameters
        ----------
        area: float
          Area of generated maps, in deg^2

        sigma: ndarray or None
          Map instrument noise, in Jy.  If None, no instrument
          noise is added.

        verbose: bool
          Print informational messages as it runs.

        Returns
        -------
          A tuple containing the input maps.  If truthtable is
        set on initialization, also includes the truth table of
        positions and fluxes, where the positions are relative
        to the first map.

        Notes
        -----
          Remember that the returned maps follow astropy.io.fits conventions,
        so the indexing into the maps is [y, x]
        """

        if area <= 0.0:
            raise ValueError("Invalid (non-positive) area")

        if sigma is None:
            int_sigma = np.zeros(self._nbands, dtype=np.float32)
        elif type(sigma) == list:
            if len(sigma) != self._nbands:
                if len(sigma) == 1:
                    int_sigma = sigma[0] * np.ones_like(self._wave)
                else:
                    raise ValueError("Number of sigmas doesn't match number"
                                     " of wavelengths")
            else:
                int_sigma = np.asarray(sigma, dtype=np.float32)
        elif type(sigma) == np.ndarray:
            if len(sigma) != self._nbands:
                if len(sigma) == 1:
                    int_sigma = sigma[0] * np.ones_like(self._wave)
                else:
                    raise ValueError("Number of sigmas doesn't match number"
                                     " of wavelengths")
            else:
                int_sigma = sigma.astype(np.float32, copy=False)
        else:
            int_sigma = float(sigma) * np.ones_like(self._wave)

        if int_sigma.min() < 0:
            raise ValueError("Invalid (negative) instrument sigma")

        # Make the non-convolved images
        # The first step is to initialize the output maps.
        # Since we do the catalog in chunks (in typical applications
        # the catalog takes more memory than the maps), we must hold
        # all the maps in memory at once
        nextent = np.empty(self._nbands, dtype=np.int32)
        truearea = np.empty(self._nbands, dtype=np.float32)
        maps = []
        for i in range(self._nbands):
            pixarea = (self._pixsize[i] / 3600.0)**2
            nextent[i] = math.ceil(math.sqrt(area / pixarea))
            truearea[i] = nextent[i]**2 * pixarea
            maps.append(np.zeros((nextent[i], nextent[i]),
                                 dtype=np.float32))

        # Figure out how many sources to make
        truearea = truearea.mean()
        nsources_base = self._npersr * (math.pi / 180.0)**2 * truearea
        nsources = np.random.poisson(lam=nsources_base)
        if verbose:
            print("True area: {:0.2f} [deg^2]".format(truearea))
            print("Number of sources to generate: {:d}".format(nsources))

        # We do this in chunks
        if self._gensize == 0:
            # One big chunk
            nchunks = 1
            chunks = np.array([nsources], dtype=np.int64)
        else:
            # Recall this is python 3 -- floating point division
            nchunks = math.ceil(nsources / self._gensize)
            print(nchunks,nsources,self._gensize)
            chunks = self._gensize * np.ones(int(nchunks), dtype=np.int64)
            chunks[-1] = nsources - (nchunks - 1) * self._gensize
            assert chunks.sum() == nsources

        # Source generation loop
        nexgen = float(nextent[0])
        if verbose:
            print("Generating sources")
        for i, nsrc in enumerate(chunks):
            if verbose and nchunks > 1:
                print("  Doing chunk {0:d} of {1:d}".format(i+1, nchunks))

            # Generate positions in base image, uniformly distributed
            # Note these are floating point
            xpos = nexgen * np.random.rand(nsrc)
            ypos = nexgen * np.random.rand(nsrc)

            # Get fluxes (in Jy)
            cat = self._gencat_stuff.generate(nsrc, wave=self._wave)
            fluxes = cat[-1].copy()

            # Set up truth table if needed
            if self._returntruth:
                truthtable = {'x': xpos, 'y': ypos,
                              'z': cat[0], 'log10M': cat[1],
                              'sb': cat[2], 'log10sSFR': cat[3],
                              'log10Lir': cat[4], 'fluxdens': fluxes}

            # Try to be good about freeing up memory
            del cat

            # Add to first map without rescaling.
            # Note this has to happen in a for loop because multiple
            # sources can go in the -same- pixel
            cmap = maps[0]
            xf = np.floor(xpos)
            yf = np.floor(ypos)
            nx, ny = cmap.shape
            np.place(xf, xf > nx-1, nx-1)
            np.place(yf, yf > ny-1, ny-1)
            for cx, cy, cf in zip(xf, yf, fluxes[:, 0]):
                cmap[int(cy), int(cx)] += cf  # Note transpose

            # Other bands, with pixel scale adjustment
            for mapidx in range(1, self._nbands):
                posrescale = self._pixsize[0] / self._pixsize[mapidx]
                xf = np.floor(posrescale * xpos)
                yf = np.floor(posrescale * ypos)
                cmap = maps[mapidx]
                nx, ny = cmap.shape
                np.place(xf, xf > nx-1, nx-1)
                np.place(yf, yf > ny-1, ny-1)
                for cx, cy, cf in zip(xf, yf, fluxes[:, mapidx]):
                    cmap[int(cy), int(cx)] += cf  # Note transpose

            if not self._returntruth:
                del fluxes, xpos, ypos, xf, yf
            else:
                del xf, yf

        # Now image details -- convolution, instrument noise
        for mapidx in range(self._nbands):
            if verbose:
                msg = "Preparing map for wavelength {0:5.1f} um "\
                      "extent: {1:d} x {2:d}"
                print(msg.format(self._wave[mapidx], maps[mapidx].shape[0],
                                 maps[mapidx].shape[1]))
                print("  Convolving")

            beam = self.get_gauss_beam(mapidx)
            maps[mapidx] = convolve(maps[mapidx], beam, boundary='wrap')

            if int_sigma[mapidx] > 0:
                if verbose:
                    msg = "  Adding instrument noise: {:0.4f} [Jy]"
                    print(msg.format(int_sigma[mapidx]))
                maps[mapidx] += np.random.normal(scale=int_sigma[mapidx],
                                                 size=maps[mapidx].shape)

            if self._returntruth:
                maps.append(truthtable)

        plt.imshow(maps[2])
        plt.title('500 micron')
        # plt.show()
        print(len(maps[-1]['fluxdens'])) # shows the truthtable
        # print(maps[2].shape)
        return maps

if __name__ == '__main__':
    gm = genmap_gauss()
    gm.generate(0.25,verbose=True)
