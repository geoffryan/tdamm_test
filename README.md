# TDAMM Models - Prototype

An experiment to implement a common interface to a variety of theoretical models electromagnetic emission models relevant to Time Domain and Multi-Messenger Astrophysics. Work began at the MOSSAIC 2023 Workshop at GWU.

## Considerations

An electromagnetic emission model is a procedure for computing the spectral flux density $F_ν$ at desired times and frequencies given a number of model parameters.  Given $F_\nu$, it is easy to compute $F_λ$, the number flux $n_ν$, integrate over a filter to generate a magnitude, etc.

### Potential Use Cases

- Basic theoretical explorations: producing some light curves or spectra locally in a notebook, varying a few parameters.
- Use by students and non-experts: simple interface, easy-to-install dependencies, obvious failure modes, good documentation.
- Parameter estimation: computing many models with varying parameters -> model evaluation should be fast & accurate
- Online Tools: integration into server-side web pages for user-drive model exploration or quick comparison to data.  Requires the calling code to *inspect* a model for details about the parameters it expects.


### Sample Architecture

A generic `TDAMMModel` class that can be instantiated with a specific model. All slow things (network communication, file I/O) should occur during object construction.

A simple `ModelParamater` class that contains information about variable parameters (name, LaTeX symbols, bounds, units).

Concrete implemenations of `TDAMMModel` contain the methods:
- `getParameters()`: Returns a list of ModelParameter objects - the parameters that can be given to `fluxDensity()`.
- `fluxDensity(t, nu, **kwargs)`: the main driver - computes the flux for the given times (t), frequencies (nu), and model parameters.  This function should *not be slow*: no file I/O or network access.  The input variables `t` and `nu` must be `astropy` Quantities with appropriate units, can be array-like, and must be the same shape or broadcastable to eachother.  `fluxDensity` returns the flux in an array the same shape as `t` and `nu` (or their broadcasted shape if `t` and `nu`'s shapes are different).


