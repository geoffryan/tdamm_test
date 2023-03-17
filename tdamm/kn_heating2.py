import numpy as np
from .tdamm_core import TDAMMModel
from astropy import constants as c
from astropy import units as u
from kilonova_heating_rate import lightcurve


def black_body_func(temp,nu):
    hnu=c.h*nu
    kT=temp*c.k_B

    x = np.atleast_1d((hnu/kT).value)
    
    one_over_expm1 = np.zeros(x.shape)
    good_idx = (x < 700)
    one_over_expm1[good_idx] = 1.0 / np.expm1(x[good_idx])

    val=2/c.c**2 * nu**2 * hnu * one_over_expm1
    return val

class KilonovaHeatingRateModel2(TDAMMModel):

    def __init__(self, mass, velocities, opacities, n, distance, z):
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
        self.zp1=z+1
        self.distance=distance
        self.mass=mass
        self.velocities=velocities
        self.opacities=opacities
        self.n=n
        #print('len(v), len(o):',len(velocities),len(opacities))
        #print('velocities',velocities)
        #print('opacities',opacities)

    
    def getParameters(self):
        return []

        
    def fluxDensity(self, times, freqs, **kwargs):
        '''
        Times, freqs assumed taken in observer frame. If array-like must be 
        broadcastable to eachother.
        Returns the spectral flux density (F_nu)
        '''

        # shift to rest-frame times and freqs
        
        ts = times / self.zp1
        nus = freqs * self.zp1

        # compute the lightcurve from the Hotokezaka heating rates

        (L, T, r) = lightcurve(ts, self.mass, self.velocities, self.opacities,
                               self.n)

        # Black-body emissivity

        Fnu = self.zp1 * black_body_func(T, nus)\
                * np.pi * (r / self.distance)**2

        # That's all folks!

        return Fnu
    

    
