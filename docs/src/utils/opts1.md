# Options and Auxiliary Functions

## Sampler Options
These options set the number of iterations and the frequency at which data is
recorderd.  This is generically used by all samplers.
```@docs
    MDOptions(; n_iters = 10^4, n_save_iters = 1)
```

## Verlet Integrator
While the main goal of this package is to sample from NVT type ensembles, as it is needed for HMC, the Verlet integrator is included should one wish to sample from the NVE ensemble:
```@docs
    Verlet
```