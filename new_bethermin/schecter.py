import numpy as np
import math

""" Schecter mass function: dN/dV dlog M = (M/Mb)^(1-alpha) Exp(-M/Mb)"""

__all__ = ["mass_schecter"]


class mass_schecter:
    """ A class for generating samples from a Schecter function.

    Because the integral of the Schecter function diverges, this requires
    a minimum and maximum mass.  This does not include the phi_b(z)
    term in the Bethermin 12 model.
    """

    def __init__(self, log10Mb=11.2, alpha=1.3, log10Mmin=9.0, log10Mmax=13.0,
                 ninterp=2000):
        """ Initializer.

        Parameters
        ----------
        log10Mb: float
          Log10 of Mb, in solar masses

        alpha: float
          Power-law slope of low mass distribution

        log10Mmin: float
          Log10 minimum mass to generate, in solar masses

        log10Mmax: float
          Log10 maximum mass to generate, in solar masses

        ninterp: float
          Number of interpolation samples to use
        """

        from scipy.interpolate import interp1d
        from scipy.integrate import trapz

        self._Mb = float(log10Mb)
        self._alpha = float(alpha)
        if self._alpha <= 1.0:
            raise ValueError("alpha should be > 1: {:f}".format(self._alpha))
        self._Mmin = float(log10Mmin)
        self._Mmax = float(log10Mmax)
        if self._Mmin == self._Mmax:
            raise ValueError("No range between min and max mass")
        if self._Mmin > self._Mmax:
            self._Mmin, self._Mmax = self._Mmax, self._Mmin
        self._ninterp = int(ninterp)
        if self._ninterp <= 0:
            raise ValueError("Ninterp must be > 0: {:d}".format(self._ninterp))

        # LogM runs from Mmax down to Mmin because the schecter function
        # drops so quickly
        logM = np.linspace(self._Mmax, self._Mmin,  self._ninterp)
        val = self(logM)

        # Needed to understand how many total samples to generate
        self._dNdV = -trapz(val, x=logM)

        # Now form inverse cumulative array we need to generate samples
        cumsum = val.cumsum()
        cumsum -= cumsum[0]  # So that 0 corresponds to the bottom
        cumsum /= cumsum[-1]  # Normalization -> 0-1 is full range
        self._interpolant = interp1d(cumsum, logM, kind='linear')

    def __call__(self, log10M):
        """ Evaluate the schechter function

        Parameters
        ----------
        log10M: float
          log10 mass to evaluate Schecter function at, in solar masses
        """

        v = 10**(log10M - self._Mb)
        return math.log(10.0) * v**(1 - self._alpha) * np.exp(-v)

    def generate(self, ngen):
        """ Generates log10 M samples from Schecter function.

        Parameters
        ----------
        ngen: int
          Number of samples to generate

        Returns
        --------
        log10M: ndarray
          Log 10 mass values, in solar masses.
        """

        return self._interpolant(np.random.rand(ngen)).astype(np.float32)

    @property
    def log10Mb(self):
        return self._Mb

    @property
    def alpha(self):
        return self._alpha

    @property
    def log10Mmin(self):
        return self._Mmin

    @property
    def log10Mmax(self):
        return self._Mmax

    @property
    def ninterp(self):
        return self._ninterp

    @property
    def totalval(self):
        return self._totval

    @property
    def dNdV(self):
        return self._dNdV

    def __str__(self):
        """ String representation"""
        outstr = "Schecter mass function with log10Mb: {0:0.2f} "\
                 "alpha: {1:0.2f} log10Mmin: {2:0.2f} log10Mmax: {3:0.2f}"
        return outstr.format(self._Mb, self._alpha, self._Mmin,
                             self._Mmax)

