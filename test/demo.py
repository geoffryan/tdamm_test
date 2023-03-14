import numpy as np
import matplotlib.pyplot as plt
from astropy import units
import tdamm

model = tdamm.AfterglowModel()

t = np.geomspace(1.0e-2, 1.0e2, 100) * units.day
nu = 1.0e14 * units.Hz

F = model(t, nu, thetaObs=10*units.deg, z=0.5, d_L=1.0e27*units.cm)

fig, ax = plt.subplots(1, 1)

ax.plot(t, F)

ax.set(xscale='log', yscale='log')

plt.show()
