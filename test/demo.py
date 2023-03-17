import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
import tdamm

afterglow_model = tdamm.AfterglowModel()
kn_model = tdamm.KN_W21("S", 1, 0.03, 0.15, 0.03, 0.15)
                        # "../kn_lanl_data/kn_sim_cube_h5")

fig, ax = plt.subplots(1, 2, figsize=(12, 4))

t = np.geomspace(1.0e-1, 1.0e3, 1000) * u.day
nu = 1.0e14 * u.Hz

F_afterglow = afterglow_model.fluxDensity(t, nu, 
        thetaObs=20*u.deg, z=0.5, d_L=1.0e27*u.cm)

F_kn = kn_model.fluxDensity(t, nu, 
        thetaObs=20*u.deg, z=0.5, d_L=1.0e27*u.cm)


ax[0].plot(t, F_afterglow.cgs, label=r'afterglow $\nu = 10^{14}$ Hz')
ax[0].plot(t, F_kn.cgs, label=r'kilonova $\nu = 10^{14}$ Hz')
ax[0].plot(t, (F_kn + F_afterglow).cgs, label='combined')

t = np.array([1.0, 10.0, 100.0, 1000.0]) * u.day
nu = np.geomspace(1.0e8, 1.0e20, 200) * u.Hz

F_afterglow = afterglow_model.fluxDensity(t[:, None], nu[None, :], 
                            thetaObs=20*u.deg, z=0.5, d_L=1.0e27*u.cm)
F_kn = kn_model.fluxDensity(t[:, None], nu[None, :], 
                            thetaObs=20*u.deg, z=0.5, d_L=1.0e27*u.cm)

ax[1].plot(nu, (F_kn + F_afterglow)[0].cgs, label=r'$t$ = 1 d')
ax[1].plot(nu, (F_kn + F_afterglow)[1].cgs, label=r'$t$ = 10 d')
ax[1].plot(nu, (F_kn + F_afterglow)[2].cgs, label=r'$t$ = 100 d')
ax[1].plot(nu, (F_kn + F_afterglow)[3].cgs, label=r'$t$ = 1000 d')

ax[0].legend()
ax[1].legend()

ax[0].set(xscale='log', yscale='log', xlabel=r'$t$ (days)',
        ylabel=r'$F_\nu$ (erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$)')
ax[1].set(xscale='log', yscale='log', xlabel=r'$\nu$ (Hz)',
        ylabel=r'$F_\nu$ (erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$)')

plt.show()
