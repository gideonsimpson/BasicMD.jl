# BasicMD


[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://gideonsimpson.github.io/BasicMD.jl/stable) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://gideonsimpson.github.io/BasicMD.jl/dev)

A collection of basic routines for Molecular Dynamics simulations written in Julia.  These include
* Euler–Maruyama (EM)
* Leimkuhler-Matthews (LM)
* Random Walk Metropolis (RWM)
* Metropolis Adjusted Langevin (MALA)
* Brünger, Brooks and Karplus (BBK)
* Grønbech-Jensen and Farago (GJF)
* ABOBA, BAOAB
* Hamiltonian/Hybrid Monte Carlo (HMC)

This is intended to be a collection for small scale problems and code
development.  Larger problems will best be handled with
[LAMMPS](https://lammps.sandia.gov/) or [OpenMM](http://openmm.org/)

# Overview
This module has been significantly reorganized to better take advantage of
multiple dispatch.  The older versions are still available in `deprecated/` as a
part of the module `OldBasicMD`, if needed.

This module can be added with the command:
```
(@v1.XYZ) pkg> add BasicMD
```

This package was previously named `JuBasicMD`, but has since been renamed to
conform with Julia package naming recommendations.

## Sampling
Each sampler is first initialized as
```
> sampler = RWM(V, β, Δt);
```
or
```
> sampler = HMC(V, gradV!, β, M, Δt, nΔt);
```
depending on the type of sampler that is chosen.  Information about what is needed to define each sampler can be obtained in the REPL by calling `?` with the sampler name (_i.e._, `?HMC`).

Once a sampler is constructed, sampling is then performed with either `sample_trajectory!` or `sample_trajectory`.  The former performs an in place transformation on the input state, while the latter records values along the trajectory:
```
sample_trajectory!(X, sampler, options=opts);
```
and
```
Xvals = sample_trajectory(X₀, sampler);
```
For Metropolis methods, the latter form also returns the running acceptance rate,
```
Xvals, avals = sample_trajectory(X₀, sampler);
```

The `sample_observables` routine is also included to allow for the computation of the observables directly.  This avoids storing an entire trajectory, which may be computationally expensive.  Given some set of observables, `f₁`, `f₂`, etc., storing these as a tuple, `obs=(f₁,f₂)`, they are computed directly with the command X   
```
obs_vals =  sample_observables(X₀ , sampler, obs)
```
This can otpionaly take the `MDOptions` in as an argument.

## MDOptions
The number of iterations performed is determined by the optional `options`
argument.  This takes as its argument a data structure which is formatted using
the `MDOptions` function:
```
opts = MDOptions(n_iters=n_iters,n_save_iters=n_save_iters)
```
`n_iters` is the number of iterations performed by the sampler.  `n_save_iters`
is the frequency with which samples are saved.    If n_save_iters=1, every
iteration is saved.  If n_save_iters=n_iters, only the final iteration is saved.


All of these methods have two versions, `Integrator` and `Integrator!`.  The `Integrator!` routine performs an in place transformation on the starting position (and momentum, where appropriate), while `Integrator` copies over the initial condition.  `Integrator` also accepts the optional argument `return_trajectory=true/false`, which will return the entire time series data.

# Constraints (Experimental)
While likely more relevant for integrating discretized SDEs than sampling, a means of constraining the flow has been introduced using the `Constraints` structure – see the examples.

# Examples
Example codes include:

* 1D Harmonic potential
* 1D Double Well Potential
* 2D Muller Potential

These examples assume that the active environment includes the `BasicMD` package.
# Caveats
* The code assumes that the state space is vector valued.  Thus, even if the problem is one dimensional, you should have initial points and functions formatted appropriately, i.e.
```
> x0 = [1.0]
```

* The mass matrix, `M`, used in the inertial Langevin integrators and Hamiltonian methods must be diagonal and provided either as a scalar (in the isotropic case) or a vector (in the anisotropic case).  This restriction is in place for performance purposes.

* BBK is currently implemented for a slightly different version of the Langevin SDE than ABOBA/BAOAB.  BBK requires inverting the mass matrix while ABOBA/BAOAB require its square root.

* GJF is implemented in (q,p) coordinates as opposed to (x,v) coordinates.  Consequently, the mass term appears slightly differently than in the literature.

Both will sample the associated Boltzmann distribution, but the SDE trajectories will differ when `M≂̸I`.

# TO DO
Before reaching a 1.0, the goals are:
* Include GHMC samplers.  
* Add to the package registry

# Acknowledgements
This work was supported in part by the US National Science Foundation Grant DMS-1818716.


# References

1. [Free Energy Computations: A Mathematical Perspective by Lelièvre, Rousset, and Stoltz](https://www.worldscientific.com/worldscibooks/10.1142/P579)
2. [Molecular Dynamics by Leimkuhler and Matthews](https://www.springer.com/gp/book/9783319163741)
3. [A simple and effective Verlet-type algorithm for simulating Langevin dynamics by Grønbech-Jesnen and Farago](https://doi.org/10.1080/00268976.2012.760055)
