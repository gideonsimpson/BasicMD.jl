# JuBasicMD
A collection of basic routines for Molecular Dynamics simulations written in Julia.  These include
* Euler–Maruyama
* Random Walk Metropolis (RWM)
* Metropolis Adjusted Langevin (MALA)
* Brünger, Brooks and Karplus (BBK)
* ABOBA, BAOAB
* Verlet

This is intended to be simple collection for small scale problems and code development.  Larger problems will best be handled with [LAMMPS](https://lammps.sandia.gov/) or [OpenMM](http://openmm.org/)

# Overview
All of these methods have two versions, `Integrator` and `Integrator!`.  The `Integrator!` routine performs an in place transformation on the starting position (and momentum, where appropriate), while `Integrator` copies over the initial condition.  `Integrator` also accepts the optional argument `return_trajectory=true/false`, which will return the entire time series data.

# Caveats
* The mass matrix, `M`, used in the inertial Langevin integrators and Hamiltonian methods must be diagonal.  This restriction is in place for performance purposes.

* BBK is currently implemented for the SDE:
$$
\begin{split}
dq& = M^{-1}p dt\\
dp & = -\nabla V(q)dt - \gamma M^{-1}p dt + \sqrt{2\beta^{-1}\gamma}dW
\end{split}
$$
In contrast, the ABOBA, BAOAB, are implemented for
$$
\begin{split}
dq& = M^{-1}p dt\\
dp & = -\nabla V(q)dt - \gamma p dt + \sqrt{2\beta^{-1}\gamma}\sqrt{M}dW
\end{split}
$$
Both will sample the associated Boltzmann distribution, but the SDE trajectories will differ when $M\neq \mathrm{Id.}$

# TO DO
Before reaching a 1.0, the goals are:
* Include HMC and GHMC samplers.
* Include reporter functions that allow for the computation of observables on at particular time intervals.  
* Migrate to a more unified calling sequence amongst the integrators.

# References

1. [Free Energy Computations: A Mathematical Perspective by Lelièvre, Rousset, and Stoltz](https://www.worldscientific.com/worldscibooks/10.1142/P579)
2. [Molecular Dynamics by Leimkuhler and Matthews](https://www.springer.com/gp/book/9783319163741)
