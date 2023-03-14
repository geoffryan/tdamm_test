import numpy as np
import matplotlib.pyplot as plt
from astropy import units
import tdamm

model = tdamm.AfterglowModel()

t = np.geomspace(1.0e-2, 1.0e2, 100) * units.day
nu = 1.0e14 * units.Hz

F = model(t, nu, jetType=0, specType=0, E0=1.0e53, n0=1.0, thetaCore=0.1, thetaWing=0.4, thetaObs=0.5, p=2.2, epsilon_e=0.1, epsilon_B=0.0001, xi_N=1.0e-2, d_L=1.0e27, z=0.5)

fig, ax = plt.subplots(1, 1)

ax.plot(t, F)

ax.set(xscale='log', yscale='log')

plt.show()
