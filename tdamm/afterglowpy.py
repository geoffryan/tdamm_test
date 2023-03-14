from .tdamm_core import TDAMMModel
import afterglowpy as grb

class AfterglowModel(TDAMMModel):

    def __init__(self):
        pass

    def fluxDensity(self, t, nu, *args, **kwargs):

        return grb.fluxDensity(t.to_value('s'), nu.to_value('Hz'),
                               *args, **kwargs)
