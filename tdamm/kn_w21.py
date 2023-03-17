from .tdamm_core import TDAMMModel, ModelParameter
from pathlib import Path
import numpy as np
import h5py as h5
import afterglowpy as grb
from astropy import units as u
from astropy import constants

class KN_W21(TDAMMModel):

    def __init__(self, topo, wind, md, vd, mw, vw, base_dir=None):
        if topo not in ["S", "P"]:
            raise ValueError("topo must be 'S' or 'P'. given: {}".format(topo))
        if wind not in [1, 2]:
            raise ValueError("wind must be 1 or 2. given: {}".format(wind))

        if base_dir is None:
            base_dir = Path(__file__).parent / "data/lanl_kn_sim"

        filename = Path(
    "LANL_KNSpec_T{0:s}_wind{1:d}_md{2:.3f}_vd{3:.3f}_mw{4:.3f}_vw{5:.3f}.h5"
        .format(topo, wind, md, vd, mw, vw))

        self._load_from_file(Path(base_dir) / filename)

    def getParameters(self):
        return [ModelParameter("thetaObs", r"$\theta_{\text{obs}}$",
                               "deg", 0.0, 180.0),
                ModelParameter("z", r"$z$",
                               "", 0.0, 20.0),
                ModelParameter("d_L", r"$d_L$",
                               "cm", 1.0e10, 1.0e30)]

    def fluxDensity(self, t_obs, nu_obs, **kwargs):

        pars = self.parseArgs(**kwargs)

        z = pars['z']
        dL = pars['d_L'].to('cm')
        thetaObs = pars['thetaObs'].to_value('rad')

        theta_idx = np.searchsorted(self.theta_bins, thetaObs) - 1
        if theta_idx < 0:
            theta_idx = 0
        if theta_idx >= self.theta_bins.shape[0]-1:
            theta_idx = self.theta_bins.shape[0]-2

        t = t_obs.to_value('s') / (1+z)
        nu = nu_obs.to_value('Hz') * (1+z)

        lam = constants.c.to_value('cm/s') / nu

        lam_idx = np.searchsorted(self.lam_bins, lam) - 1

        t_idx = np.searchsorted(self.t, t) - 1

        b = np.broadcast(t, t_idx, lam_idx)

        Nt = self.t.shape[0]
        Nlam = self.lam_bins.shape[0]-1

        Llam = np.empty(b.shape)
        Llam.flat = [((self.t[i+1]-tval) * self.L_lam[i-1, j, theta_idx]
                      + (tval - self.t[i]) * self.L_lam[i, j, theta_idx]
                     ) / (self.t[i+1] - self.t[i])
                        if (i >= 0 and i < Nt-1 and j >= 0 and j < Nlam)
                        else 0.0
                    for (tval, i, j) in b]


        Lnu = (constants.c / (nu * u.Hz)**2) * (Llam * u.erg/(u.s * u.cm))

        return ((1+z) * Lnu / (4*np.pi * dL**2)).to('erg/(s cm2 Hz)')


    def _load_from_file(self, path):

        with h5.File(path) as f:
            self.topo = "P" if f['topo'][...] == 1 else "S"
            self.wind = f['wind'][...]
            
            self.md = f['md_Msolar'][...]
            self.vw = f['vw_c'][...]
            self.mw = f['mw_Msolar'][...]
            self.vd = f['vd_c'][...]
            
            self.t = (f['t_days'][...] * u.day).to_value('s')
            lam_e = f['lambda_cm'][...]
            theta_e = f['theta_rad'][...]
            self.L_lam = (f['fla_cgs_per_angstrom'][...]
                          * u.erg / (u.s * u.cm**2 * u.Angstrom)
                          * 54 * 4*np.pi * (10.0 * u.pc)**2).to_value(
                                  'erg/(s cm)')

        self.lam_bins = np.empty((lam_e.shape[0]+1), dtype=lam_e.dtype)
        self.theta_bins = np.empty((theta_e.shape[0]+1), dtype=theta_e.dtype)

        self.lam_bins[:-1] = lam_e[:, 0]
        self.lam_bins[-1] = lam_e[-1, 1]
        self.theta_bins[:-1] = theta_e[:, 0]
        self.theta_bins[-1] = theta_e[-1, 1]
