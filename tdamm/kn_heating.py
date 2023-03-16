
from .tdamm_core import TDAMMModel
import afterglowpy as grb


#!/usr/bin/env python
import timeit

from astropy import constants as c
from astropy import units as u
from kilonova_heating_rate import lightcurve


mass = 0.05 * u.Msun
velocities = np.asarray([0.1, 0.2, 0.4]) * c.c
opacities = np.asarray([3.0, 0.5]) * u.cm**2 / u.g
n = 4.5
t = np.geomspace(0.02, 10) * u.day

L, T, r = lightcurve(t, mass, velocities, opacities, n)

# Evaluate flux in band at 100 Mpc for a few different filters.
DL = 100 * u.Mpc
seds = [
    synphot.SourceSpectrum(synphot.BlackBody1D, temperature=TT)
    * np.pi * (rr / DL).to(u.dimensionless_unscaled)**2
    for TT, rr in zip(T, r)]

def black_body_func(temp,nu):
    hnu=c.h*nu
    kT=temp*c.k
    val=2/c.c**2*nu**2*hnu/np.expm1(hnu/kT)
    return val

class KilonovaModel(TDAMMModel):

    def __init__(self, mass, velocities, opacities, n, distance, z):
        'Kilonova heating model
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
        
        Ref
        '''
        self.zp1=z+1
        self.distance=distance
        self.mass=mass
        self.velocities=velocities
        self.opacities=opacities
        self.n=n


    def TFgrids(self, times, freqs):
        '''
        Identify gridded data. This should go into the base class
        '''

        if times.shape[0]==1 and freqs.shape[0]==1:
            return (times[0],freqs[0])
        else: return None        

        
    def fluxDensity(self, times, freqs, *args, **kwargs):
        '''
        Times, freqs assumed taken in observer frame. The challenge is whether        
        
        '''
        
        nus=freqs*zp1

        tf = self.TFgrids(times,freqs)
        if tf is None:            
            # Evaluate flux in band at 100 Mpc for a few different filters.
            ts=times/self.zp1
            nus=freqs*self.zp1
            (L, T, r) = lightcurve(ts, self.mass,self.velocities, self.opacities, self.n)            
            seds=[ black_body_func(TT,nu) * np.pi * self.zp1 (rr / self.distance).to(u.dimensionless_unscaled)**2
                   for TT, rr, nu in zip(T, r, nus) ]
            return np.array(seds)
        else:
            ts=times[0]/self.zp1
            nus=freqs[0]*self.zp1
            (L, T, r) = lightcurve(ts, self.mass,self.velocities, self.opacities, self.n)            
            seds=[ black_body_func(TT,nus) * np.pi * self.zp1 (rr / self.distance).to(u.dimensionless_unscaled)**2
                   for TT, rr in zip(T, r) ]
            return np.array(seds)[None]

    

    
