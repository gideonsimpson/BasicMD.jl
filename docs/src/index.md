# BasicMD.jl Documentation


## Overview
This is a collection of basic routines for Molecular Dynamics simulations
written in Julia.  These include
* Euler–Maruyama (EM)
* Leimkuhler-Matthews (LM)
* Random Walk Metropolis (RWM)
* Metropolis Adjusted Langevin (MALA)
* Brünger, Brooks and Karplus (BBK)
* Grønbech-Jensen and Farago (GJF)
* ABOBA, BAOAB
* Hamiltonian/Hybrid Monte Carlo (HMC)

## Caveats
* The code assumes that the state space is vector valued.  Thus, even if the
  problem is one dimensional, you should have initial points and functions
  formatted appropriately, i.e.
```
> x0 = [1.0]
```

* The mass matrix, `M`, used in the inertial Langevin integrators and
  Hamiltonian methods must be diagonal and provided either as a scalar (in the
  isotropic case) or a vector (in the anisotropic case).  This restriction is in
  place for performance purposes.

* BBK is currently implemented for a slightly different version of the Langevin
  SDE than ABOBA/BAOAB.  BBK requires inverting the mass matrix while
  ABOBA/BAOAB require its square root.

* GJF is implemented in (q,p) coordinates as opposed to (x,v) coordinates.
  Consequently, the mass term appears slightly differently than in the
  literature.

Both will sample the associated Boltzmann distribution, but the SDE trajectories
will differ when `M≂̸I`.

## Acknowledgements
This work was supported in part by the US National Science Foundation Grant DMS-1818716.
