from .tdamm_core import TDAMMModel, ModelParameter
import afterglowpy as grb

class AfterglowModel(TDAMMModel):

    def __init__(self):
        pass

    def getParameters(self):
        return [ModelParameter("thetaObs", r"$\theta_{\text{obs}}$",
                               "deg", 0.0, 180.0),
                ModelParameter("z", r"$z$",
                               "", 0.0, 20.0),
                ModelParameter("d_L", r"$d_L$",
                               "cm", 1.0e10, 1.0e30)]

    def fluxDensity(self, t, nu, **kwargs):

        pars = self.parseArgs(**kwargs)

        Z = dict(jetType=grb.jet.Gaussian, specType=0,
                 E0=1.0e53, n0=1.0e-3, thetaCore=0.06, thetaWing=0.4,
                 p=2.2, epsilon_e=0.1, epsilon_B=1.0e-3, xi_N=1.0e-2,
                 thetaObs=pars['thetaObs'].to_value("rad"), 
                 z=pars['z'],
                 d_L=pars['d_L'].to_value('cm'))

        return grb.fluxDensity(t.to_value('s'), nu.to_value('Hz'),
                               **Z)
