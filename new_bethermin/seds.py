import numpy as np
import math

""" Handles SED templates for Bethermin et al. 2012 model"""

__all__ = ["sed_model"]


class sed_model:
    """ Gets SEDs and fluxes from Magdis template models"""

    def __init__(self, Om0=0.315, H0=67.7, zmin=0.5,
                 zmax=7.0, ninterp=100):
        """ Initializer

        Parameters
        ----------
        Om0: float
          Matter density parameter

        H0: float
          Hubble constant, in km / s / Mpc

        zmin: float
          Minimum redshift supported.

        zmax: float
          Maximum redshift supported.

        ninterp: int
          Number of interpolation points for luminosity distance.

        """

        from astropy.cosmology import FlatLambdaCDM
        from astropy.units import Quantity
        import astropy.io.fits as fits
        from pkg_resources import resource_filename
        from scipy.interpolate import RectBivariateSpline as rbs
        from scipy.interpolate import interp1d

        self._zmin = float(zmin)
        self._zmax = float(zmax)
        self._Om0 = float(Om0)
        self._H0 = float(H0)
        self._ninterp = int(ninterp)

        if self._zmin == self._zmax:
            raise ValueError("No range between zmin and zmax")
        if self._zmin > self._zmax:
            self._zmin, self._zmax = self._zmax, self._zmin
        if self._zmin < 0.0:
            raise ValueError("zmin must be >= 0: {:f}".format(self._zmin))
        if self._Om0 <= 0.0:
            raise ValueError("Om0 must be positive: {:f}".format(self._Om0))
        if self._H0 <= 0.0:
            raise ValueError("H0 must be positive: {:f}".format(self._H0))
        if self._ninterp <= 0:
            raise ValueError("Ninterp must be > 0: {:d}".format(self._ninterp))

        # Set up luminosity distance interpolant.  We actually
        # interpolate log((1+z) / (4 pi d_L^2)) in log(1+z)
        cos = FlatLambdaCDM(H0=self._H0, Om0=self._Om0, Tcmb0=0.0, Neff=0.0)
        zrange = np.linspace(self._zmin, self._zmax, self._ninterp)
        mpc_in_cm = 3.0857e24
        prefac = 1.0 / (4 * math.pi * mpc_in_cm**2)
        lumdist = cos.luminosity_distance(zrange)
        if isinstance(lumdist, Quantity):
            lumdist = lumdist.value
        dlval = prefac * (1.0 + zrange) / lumdist**2
        self._dlfac = interp1d(np.log(1 + zrange), np.log(dlval))

        # Read in the data products, and set up interpolations on them
        sb_tpl = resource_filename(__name__, 'resources/SED_sb.fits')
        hdu = fits.open(sb_tpl)
        dat = hdu['SEDS'].data
        hdu.close()
        self._sblam = dat['LAMBDA'][0]
        self._sbumean = dat['UMEAN'][0]
        arg = np.argsort(self._sbumean)
        self._sbumean = self._sbumean[arg]
        self._sbrange = np.array([self._sbumean[0], self._sbumean[-1]])
        self._sbseds = dat['SEDS'][0, :, :].transpose()[arg, :]
        self._sbinterp = rbs(self._sbumean, self._sblam, self._sbseds,
                             kx=1, ky=1)

        ms_tpl = resource_filename(__name__, 'resources/SED_ms.fits')
        hdu = fits.open(ms_tpl)
        dat = hdu['SEDS'].data
        hdu.close()
        self._mslam = dat['LAMBDA'][0]
        self._msumean = dat['UMEAN'][0]
        arg = np.argsort(self._msumean)
        self._msumean = self._msumean[arg]
        self._msrange = np.array([self._msumean[0], self._msumean[-1]])
        self._msseds = dat['SEDS'][0, :, :].transpose()[arg, :]
        self._msinterp = []
        self._msinterp = rbs(self._msumean, self._mslam, self._msseds,
                             kx=1, ky=1)

    @property
    def zmin(self):
        return self._zmin

    @property
    def zmax(self):
        return self._zmax

    def get_sed(self, z, U, is_starburst, log10lir=1.0):
        """ Gets the combined SED in erg/s/cm^2/Hz"""

        zval = float(z)
        if zval > self._zmax or zval < self._zmin:
            raise ValueError("z out of supported range: {:f}".format(zval))
        opz = 1.0 + zval
        ldfac = 10**log10lir * np.exp(self._dlfac(np.log(opz)))

        if is_starburst:
            return (self._sblam / opz,
                    ldfac * self._intsed1(U, self._sbumean, self._sbseds))
        else:
            return (self._mslam / opz,
                    ldfac * self._intsed1(U, self._msumean, self._msseds))

    def _intsed1(self, U, uarr, seds):
        # uarr[idx] <= U < uarr[idx+1]
        idx = np.searchsorted(uarr, U, side='right')
        if idx == 0:
            return seds[0, :]
        elif idx == len(uarr):
            return seds[-1, :]
        wt2 = (U - uarr[idx-1])/(uarr[idx] - uarr[idx-1])
        wt1 = 1.0 - wt2
        return wt1 * seds[idx-1, :] + wt2 * seds[idx, :]

    def get_fluxes(self, wave, z, U, is_starburst, log10lir=0.0):
        """ Gets the flux density at the observer frame wavelengths wave
        in Jy for a set of sources with a given log10 L_IR.  They must
        be all starburst, or all main sequence"""

        # Simple scalar case
        if np.isscalar(z):
            if not np.isscalar(U):
                raise ValueError("U must be scalar if z is")
            if not np.isscalar(log10lir):
                raise ValueError("log10lir must be scalar if z is")
            zval = float(z)
            if zval > self._zmax:
                errmsg = "z {0:f} above supported range: {1:f}"
                raise ValueError(errmsg.format(zval, self._zmax))
            if zval < self._zmin:
                errmsg = "z {0:f} below supported range: {1:f}"
                raise ValueError(errmsg.format(zval, self._zmin))
            opz = 1.0 + zval
            ldfac = 10**(log10lir + 23.0) * \
                np.exp(self._dlfac(np.log(opz)))  # 1e23 is to Jy

            if is_starburst:
                rng = self._sbrange
                interp = self._sbinterp
            else:
                rng = self._msrange
                interp = self._msinterp
            if U < rng[0]:
                U = rng[0]
            elif U > rng[1]:
                U = rng[1]
            if np.isscalar(wave):
                return ldfac * interp(U, wave / opz).flatten()
            else:
                return ldfac * interp(U, np.asarray(wave) / opz).flatten()
        else:
            # Not so simple case
            z = np.asarray(z)
            nz = len(z)
            if np.isscalar(U):
                U = U * np.ones_like(z)
            else:
                U = np.asarray(U)
                if len(U) != nz:
                    raise ValueError("Mismatch between number of z and "
                                     "U values")
            if np.isscalar(log10lir):
                log10lir = log10lir * np.ones_like(z)
            else:
                log10lir = np.asarray(log10lir)
                if len(log10lir) != nz:
                    raise ValueError("Mismatch between number of z and "
                                     "log10lir values")

            if z.max() > self._zmax:
                errmsg = "z {0:f} above supported range: {1:f}"
                raise ValueError(errmsg.format(z.max(), self._zmax))
            if z.min() < self._zmin:
                errmsg = "z {0:f} below supported range: {1:f}"
                raise ValueError(errmsg.format(z.min(), self._zmin))
            opz = 1.0 + z
            ldfac = 10**(log10lir + 23.0) * \
                np.exp(self._dlfac(np.log(opz)))  # 1e23 is to Jy

            if is_starburst:
                rng = self._sbrange
                interp = self._sbinterp
            else:
                rng = self._msrange
                interp = self._msinterp
            np.copyto(U, rng[0], where=(U < rng[0]))
            np.copyto(U, rng[1], where=(U > rng[1]))

            if np.isscalar(wave):
                retarr = np.empty_like(z)
                for idx in range(nz):
                    retarr[idx] = ldfac[idx] *\
                        interp(U[idx], wave / opz[idx]).flatten()
            else:
                retarr = np.empty((nz, len(wave)), dtype=np.float32)
                for idx in range(nz):
                    retarr[idx, :] = ldfac[idx] *\
                        interp(U[idx], wave / opz[idx]).flatten()
            return retarr
