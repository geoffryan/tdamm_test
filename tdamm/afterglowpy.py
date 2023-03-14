from .tdamm_core import TDAMMModel, ModelParameter
import afterglowpy as grb

class AfterglowModel(TDAMMModel):

    def __init__(self):
        pass

    def getParameters(self):
        return [ModelParameter("thetaObs", r"$\theta_{\text{obs}}$",
                               "deg", 0.0, 180.0)]

    def fluxDensity(self, t, nu, *args, **kwargs):

        return grb.fluxDensity(t.to_value('s'), nu.to_value('Hz'),
                               *args, **kwargs)
