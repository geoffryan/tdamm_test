import numpy as np
from astropy import units

class TDAMMModel:

    def __init__(self):
        pass

    def fluxDensity(self, t, nu, **kwargs):
        pass

    def getParameters(self):
        pass

    def __call__(self, t, nu, **kwargs):
        return self.fluxDensity(t, nu, **kwargs)

    def parseArgs(self, **kwargs):

        pars = self.getParameters()

        Z = {}
        for par in pars:
            if par.name not in kwargs:
                raise RuntimeError("Parameter {0:s} not given to model"
                                    .format(par.name))
            Z[par.name] = kwargs[par.name]

        return Z


class ModelParameter:

    def __init__(self, name, name_tex, unit_name, low, high):
        self.name = name
        self.name_tex = name_tex
        self.unit_name = unit_name
        self.low = low
        self.high = high


class Dummy(TDAMMModel):

    def __init__(self):
        pass

    def fluxDensity(self, t, nu, F0, alpha, beta):
        t1d = 1.0 * units.day
        nu0 = 1.0e14 * units.Hz

        return F0 * np.power(t/t1d, alpha) * np.power(nu/nu0, beta)
