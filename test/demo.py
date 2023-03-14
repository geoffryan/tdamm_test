import numpy as np
import matplotlib.pyplot as plt
from astropy import units
import tdamm

model = tdamm.Dummy()

t = np.geomspace(1.0e-2, 1.0e2, 100) * units.day
nu = np.geomspace(1.0e6, 1.0e14, 100) * units.Hz

F = model(t, nu, 1.0, 0.5, -2.5)

fig, ax = plt.subplots(1, 1)

ax.plot(nu, F)

ax.set(xscale='log', yscale='log')

plt.show()
