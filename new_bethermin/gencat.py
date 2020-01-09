import numpy as np
import math

""" Generates sources from the Bethermin et al. 2012 model"""

# __all__ = ["gencat"]

from schecter import mass_schecter
import zdist
from seds import sed_model


class Gencat:
    """ Generates catalog sources from Bethermin et al. 2012 model"""

    def __init__(self, log10Mb=11.2, alpha=1.3, log10Mmin=8.0, log10Mmax=12.75,
                 ninterpm=2000, zmin=0.1, zmax=10.0, Om0=0.315,
                 H0=67.7, phib0=-3.02, gamma_sfmf=0.4, ninterpz=1000,
                 rsb0=0.012, gammasb=1.0, zsb=1.0, logsSFRM0=-10.2,
                 betaMS=-0.2, zevo=2.5, gammams=3.0, bsb=0.6, sigmams=0.15,
                 sigmasb=0.2, mnU_MS0=4.0, gammaU_MS0=1.3, z_UMS=2.0,
                 mnU_SB0=35.0, gammaU_SB0=0.4, z_USB=3.1, scatU=0.2,
                 ninterpdl=200):
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

        self._zmin = float(zmin)
        self._zmax = float(zmax)
        self._rsb0 = float(rsb0)
        self._gammasb = float(gammasb)
        self._zsb = float(zsb)
        self._interp = int(ninterpz)
        self._logsSFRM0 = float(logsSFRM0)
        self._betaMS = float(betaMS)
        self._zevo = float(zevo)
        self._gammams = float(gammams)
        self._bsb = float(bsb)
        self._sigmams = float(sigmams)
        self._sigmasb = float(sigmasb)
        self._scatU = float(scatU)  # dex
        self._scatUe = math.log(10.0) * self._scatU  # log_10 -> ln
        self._mnUMS = float(mnU_MS0)
        self._gammaUMS = float(gammaU_MS0)
        self._zUMS = float(z_UMS)
        self._mnUSB = float(mnU_SB0)
        self._gammaUSB = float(gammaU_SB0)
        self._zUSB = float(z_USB)
        self._Om0 = float(Om0)
        self._H0 = float(H0)

        if self._mnUMS <= 0:
            errstr = "mnU_MS0 must be positive, not {:f}".format(self._mnUMS)
            raise ValueError(errstr)
        if self._zUMS < 0:
            errstr = "z_UMS must be non-negative, not {:f}".format(self._zUMS)
            raise ValueError(errstr)
        if self._mnUSB <= 0:
            errstr = "mnU_SB0 must be positive, not {:f}".format(self._mnUSB)
            raise ValueError()
        if self._zUSB < 0:
            errstr = "z_USB must be non-negative, not {:f}".format(self._zUSB)
            raise ValueError(errstr)

        self._sch = mass_schecter(log10Mb, alpha, log10Mmin, log10Mmax,
                                  ninterpm)
        self._zdist = zdist.Zdist(self._zmin, self._zmax, self._Om0, self._H0,
                            phib0, gamma_sfmf, ninterpz)
        self._sedmodel = sed_model(zmin=self._zmin, zmax=self._zmax,
                                   Om0=self._Om0, H0=self._H0,
                                   ninterp=ninterpdl)

        # Set number per sr
        self._npersr = self._zdist.dVPhidzdOmega * self._sch.dNdV

    @property
    def npersr(self):
        return self._npersr

    @property
    def npersqdeg(self):
        return self._npersr * (math.pi / 180.0)**2

    def generate(self, ngen, wave=None):
        """ Generates samples from the Bethermin 2012 model.

        Returns a tuple of (z, log10 M, is_starburst, log10 sSFR),
        each of which is a ngen element ndarray.  If wave is
        not None, will also generate flux densities (in Jy) for
        each source."""

        log10mass = self._sch.generate(ngen)
        z = self._zdist.generate(ngen)

        # Figure out if each source is a starburst
        rsb = (1.0 + self._zsb)**self._gammasb * np.ones(ngen)
        w = np.nonzero(z < self._zsb)[0]
        if len(w) > 0:
            rsb[w] = (1.0 + z[w])**self._gammasb
        rsb *= self._rsb0
        prob_starburst = self._sigmasb * rsb /\
                         (self._sigmams + self._sigmasb * rsb)
        is_starburst = np.random.rand(ngen) < prob_starburst
        del rsb, prob_starburst

        # Figure out sSFR for each source.  These are gaussian -- just
        # times different numbers and means depending on whether they
        # are a SB, plus redshift, plus mass
        # Redshift evolution of MS value

        logsSFRM = self._logsSFRM0 + self._betaMS * (log10mass - 11.0)
        w = np.nonzero(z < self._zevo)[0]
        if len(w) > 0:
            logsSFRM[w] += self._gammams * np.log10(1.0 + z[w])
        w = np.nonzero(z >= self._zevo)[0]
        if len(w) > 0:
            logsSFRM[w] += self._gammams * math.log10(1.0 + self._zevo)

        # Setup sigmas as well
        sigmas = self._sigmams * np.ones(ngen)
        w = np.nonzero(is_starburst)[0]
        if len(w) > 0:
            logsSFRM[w] += self._bsb
            sigmas[w] = self._sigmasb

        # Actual generation of log sSFR
        log10sSFR = logsSFRM + sigmas * np.random.randn(ngen)
        del logsSFRM
        del sigmas
        del w

        print('wave:',wave)
        if not wave is None:
            nwave = len(wave)
            fluxes = np.empty((ngen, nwave), dtype=np.float32)
            kfac = math.log10(1.7e-10)  # Kennicutt '98 conversion

            # Get log10 lir; note I'm assuming the r1500 business
            # only applies to SBs, not MSs
            log10lir = (log10mass + log10sSFR - kfac).astype(np.float32)

            # Deal with extinction effects on L_IR
            # coeff values are from eq 7 of B12 * 0.4 (from eq 8)
            pow_r1500 = 10**(1.628 * log10mass - 15.728)
            fsf = pow_r1500 / (1.0 + pow_r1500)
            log10lir += np.log10(fsf).astype(np.float32)

            # Do starbursts
            wsb = np.nonzero(is_starburst)[0]
            nsb = len(wsb)
            if nsb > 0:
                # Get the U (mean radiation field), eq 6 in B12
                u = self._mnUSB * \
                    (1.0 + np.minimum(z[wsb], self._zUSB))**self._gammaUSB
                # Add scatter to U
                if self._scatU > 0.0:
                    u *= np.random.lognormal(sigma=self._scatUe, size=(nsb))
                # Actual fluxes
                fluxes[wsb, :] = \
                     self._sedmodel.get_fluxes(wave, z[wsb], u, True,
                                               log10lir=log10lir[wsb])
            del wsb

            # Do MS
            wms = np.nonzero(~is_starburst)[0]
            nms = len(wms)
            if nms > 0:
                u = self._mnUMS * \
                     (1.0 + np.minimum(z[wms], self._zUMS))**self._gammaUMS
                if self._scatU > 0.0:
                    u *= np.random.lognormal(sigma=self._scatUe, size=(nms))
                fluxes[wms, :] = \
                      self._sedmodel.get_fluxes(wave, z[wms], u, False,
                                                log10lir=log10lir[wms])
            del wms

            return (z, log10mass, is_starburst, log10sSFR, log10lir, fluxes)
        else:
            return (z, log10mass, is_starburst, log10sSFR)
