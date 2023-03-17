import numpy as np
from .tdamm_core import TDAMMModel, ModelParameter
from astropy import constants as c
from astropy import units as u
from kilonova_heating_rate import lightcurve


def black_body_func(temp,nu):
    # The Planck function
    
    hnu=c.h*nu
    kT=temp*c.k_B

    # This hacky thing is to stop overflow warnings.
    x = np.atleast_1d((hnu/kT).value)    # put argument in a clean numpy array
    one_over_expm1 = np.zeros(x.shape)   # initialize destination arr to 0.
    good_idx = (x < 700)                 # find the values that won't overflow
    one_over_expm1[good_idx] = 1.0 / np.expm1(x[good_idx])  # compute those!

    # now "one_over_expm1" contains the correct values of 1/(e^x - 1).
    # The skipped values, x > 700 ==> e^x > 10^304  ==> 1/(e^x-1) < 10^(-304)
    # that is, we truncated function values below 10^(-304) to 0.0, which is
    # fine.

    # the actual planck function now!
    val = 2/c.c**2 * nu**2 * hnu * one_over_expm1

    return val

class KilonovaHeatingRateModel2(TDAMMModel):

    def __init__(self, mass, velocities, opacities, n):
        '''
        Kilonova heating rate model

        args:
          velocities : :class:`astropy.units.Quantity`
            Array of ejecta velocities in units compatible with `cm/s`.
            Length must be >= 2.
          opacities : :class:`astropy.units.Quantity`
            Array of opacities in units compatible with `cm**2/g`.
            Lenght must be >= 1, and 1 less than the length of `vej`.
          n : int, float
            Power-law index of density profile.
          distance: float
            Luminosity distance to source
          z: float
            redshift
        
        Ref:
        https://github.com/Basdorsman/kilonova-heating-rate

        If you use this work to produce a peer-reviewed journal article, please cite the following papers:

        Korobkin, O., Rosswog, S., Arcones, A., & Winteler, C. 2012, "On the astrophysical robustness of the neutron star merger r-process," Monthly Notices of the Royal Astronomical Society, 426, 1940. https://doi.org/10.1111/j.1365-2966.2012.21859.x
        Hotokezaka, K. & Nakar, E. 2020, "Radioactive Heating Rate of r-process Elements and Macronova Light Curve," Astrophysical Journal, 891, 152. https://doi.org/10.3847/1538-4357/ab6a98
        '''
        self.mass=mass
        self.velocities=velocities
        self.opacities=opacities
        self.n=n
        #print('len(v), len(o):',len(velocities),len(opacities))
        #print('velocities',velocities)
        #print('opacities',opacities)

    
    def getParameters(self):
        return [ModelParameter("z", r"$z$",
                               "", 0.0, 20.0),
                ModelParameter("d_L", r"$d_L$",
                               "cm", 1.0e10, 1.0e30)]

        
    def fluxDensity(self, times, freqs, **kwargs):
        '''
        Times, freqs assumed taken in observer frame. If array-like must be 
        broadcastable to eachother.
        Returns the spectral flux density (F_nu)
        '''

        # First, we grab our runtime arguments from the kwargs.
        # We do this like this so we only have to put the argument list in
        # one place (getParameters())
        pars = self.parseArgs(**kwargs)

        # "pars" is now a dict containing our arguments, how nice!
        z = pars['z']
        dL = pars['d_L'].to('cm')
        

        # shift to rest-frame times and freq
        ts = times / (1 + z)
        nus = freqs * (1 + z)

        # compute the temperatue and radius from the Hotokezaka heating rates
        (L, T, r) = lightcurve(ts, self.mass, self.velocities, self.opacities,
                               self.n)

        # Black-body emissivity
        Fnu = (1 + z) * black_body_func(T, nus) * np.pi * (r / dL)**2

        # That's all folks!

        return Fnu
    

    
