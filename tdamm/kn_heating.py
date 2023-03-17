import numpy as np
from .tdamm_core import TDAMMModel
#import synphot

#!/usr/bin/env python

from astropy import constants as c
from astropy import units as u
from kilonova_heating_rate import lightcurve


def black_body_func(temp,nu):
    hnu=c.h*nu
    kT=temp*c.k_B
    val=2/c.c**2*nu**2*hnu/np.expm1((hnu/kT).to(u.dimensionless_unscaled))
    return val

class KilonovaHeatingRateModel(TDAMMModel):

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


    def TFgrids(self, times, freqs):
        '''
        Identify gridded data. This should go into the base class
        '''

        if times.shape[0]==1 and freqs.shape[0]==1:
            return (times[0],freqs[0])
        else: return None        

        
    def fluxDensity(self, times, freqs, *args, **kwargs):
        '''
        Times, freqs assumed taken in observer frame. The challenge is whether. 
        Any times < 0.02 day are ignored.
        
        '''
        times=np.atleast_1d(times)
        tf = self.TFgrids(times,freqs)
        if tf is None:            
            # Evaluate flux in band at 100 Mpc for a few different filters.
            ts=times/self.zp1
            #keep=ts>0.02*u.d
            #ts=ts[keep]
            nus=freqs*self.zp1
            if len(nus.shape)==0: nus=nus*np.ones(len(ts))
            #print('ts',ts)
            #print('nus',nus)
            (L, T, r) = lightcurve(ts, self.mass, self.velocities, self.opacities, self.n)
            sedunits=u.erg/u.cm**2
            seds=[ (
                    ( black_body_func(TT,nu) * np.pi * self.zp1* (rr / self.distance)**2 ) /sedunits
                   ).to(u.dimensionless_unscaled)
                for TT, rr, nu in zip(T, r, nus) ]
            seds= np.asarray(seds)*sedunits
            
        else:
            ts=times[0]/self.zp1
            #keep=ts>0.02*u.d
            #ts=ts[keep]
            nus=freqs[0]*self.zp1
            (L, T, r) = lightcurve(ts, self.mass, self.velocities, self.opacities, self.n)            
            seds=[ black_body_func(TT,nus) * np.pi * self.zp1 * (rr / self.distance).to(u.dimensionless_unscaled)**2
                   for TT, rr in zip(T, r) ]
            seds=np.array(seds)[None]
        return seds
    

    
