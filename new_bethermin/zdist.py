import numpy as np
import math

""" Redshift distribution for Bethermin et al. 2012 model.

This contains the comoving volume element and the phi_b(z) function.
"""

__all__ = ["zdist"]


class zdist:
    """ Volume element and z distribution"""

    def __init__(self, zmin=0.5, zmax=7.0, Om0=0.315, H0=67.7, phib0=-3.02,
                 gamma_sfmf=0.4, ninterp=1000):
        """ Initializer.

        Parameters
        ----------
        zmin : float
          Minimum redshift

        zmax: float
          Maximum redshift

        Om0: float
          Matter density parameter

        H0: float
          Hubble constant in km / s / Mpc

        phib0: float
          log 10 number density at SFMF break in comoving Mpc^-3

        gamma_sfmf: float
          Evolution of the number density of the SFMF at z > 1

        ninterp: int
          Number of interpolation samples to use

        """

        from scipy.interpolate import interp1d
        from scipy.integrate import trapz
        from astropy.cosmology import FlatLambdaCDM
        from astropy.units import Quantity

        self._zmin = float(zmin)
        self._zmax = float(zmax)
        if self._zmin == self._zmax:
            raise ValueError("No range between zmin and zmax")
        if self._zmin > self._zmax:
            self._zmin, self._zmax = self._zmax, self._zmin
        if self._zmin < 0.0:
            raise ValueError("zmin must be >= 0: {:f}".format(self._zmin))
        self._Om0 = float(Om0)
        if self._Om0 <= 0.0:
            raise ValueError("Om0 must be positive: {:f}".format(self._Om0))
        self._H0 = float(H0)
        if self._H0 <= 0.0:
            raise ValueError("H0 must be positive: {:f}".format(self._H0))
        self._ninterp = int(ninterp)
        if self._ninterp <= 0:
            raise ValueError("Ninterp must be > 0: {:d}".format(self._ninterp))
        self._phib0 = float(phib0)
        self._gamma_sfmf = float(gamma_sfmf)

        self._zvals = np.linspace(self._zmin, self._zmax, self._ninterp)

        # Cosmology bit
        c_over_H0 = 299792.458 / self._H0  # in Mpc
        cos = FlatLambdaCDM(H0=self._H0, Om0=self._Om0)

        # in comoving Mpc^3
        ang_diam_dist = cos.angular_diameter_distance(self._zvals)
        if isinstance(ang_diam_dist, Quantity):
            ang_diam_dist = ang_diam_dist.value
        self._dVdzdOmega = c_over_H0 * (1.0 + self._zvals)**2 * \
            ang_diam_dist**2 / np.sqrt((1.0 + self._zvals)**3 *
                                       self._Om0 + (1.0 - self._Om0))

        # Schecter evolution bit
        phi = self._phib0 * np.ones(self._ninterp, dtype=np.float64)
        wgt1 = np.nonzero(self._zvals > 1.0)[0]
        if len(wgt1) > 0:
            phi[wgt1] += self._gamma_sfmf * (1.0 - self._zvals[wgt1])

        # Combined
        comb = 10**phi * self._dVdzdOmega

        # Needed to understand normalization
        self._dVPhidzdOmega = trapz(comb, x=self._zvals)

        # Form inverse cumulative array needed to generate samples
        cumsum = comb.cumsum()
        cumsum -= cumsum[0]  # So that 0 corresponds to the bottom
        cumsum /= cumsum[-1]  # Normalization -> 0-1 is full range
        self._interpolant = interp1d(cumsum, self._zvals, kind='linear')

    def generate(self, ngen):
        """ Generates z samples from redshift distribution.

        Parameters
        ----------
        ngen: int
          Number of samples to generate.

        Returns
        -------
        z: ndarray
          Redshifts
        """

        return self._interpolant(np.random.rand(ngen)).astype(np.float32)

    @property
    def zmin(self):
        return self._zmin

    @property
    def zmax(self):
        return self._zmax

    @property
    def zvals(self):
        """ Tabulated redshift values"""
        return self._zvals

    @property
    def Om0(self):
        return self._Om0

    @property
    def H0(self):
        return self._H0

    @property
    def phib0(self):
        return self._phib0

    @property
    def gamma_sfmf(self):
        return self._gamma_sfmf

    @property
    def dVdzdOmega(self):
        """ Comoving volume in Mpc^3 per redshift per sr at zvals"""
        return self._dVdzdOmega

    @property
    def dVPhidzdOmega(self):
        """ Comoving volume in Mpc^3 per redshift per sr at zvals * phi_b(z)"""
        return self._dVPhidzdOmega

