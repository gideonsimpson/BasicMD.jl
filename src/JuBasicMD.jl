"""
# JuBasicMD.jl

A collection of basic routines for Molecular Dynamics simulations written in
Julia.  These include:
* Euler–Maruyama (EM)
* Random Walk Metropolis (RWM)
* Metropolis Adjusted Langevin (MALA)
* Brünger, Brooks and Karplus (BBK)
* ABOBA, BAOAB
* Verlet
* Hamiltonian/Hybrid Monte Carlo (HMC)

## REPL help
`?` followed by an sampler name (`?RWM`)to obtain infromation about individual
samplers.  Use `?` on `sample_trajectory!` and `sample_trajectory` to obtain
information about running the samplers.  Finally, `?` on `MDOptions` will provide
information about how to set additional options, such as the number of
information.

## Performing Sampling

Having constructed a `sampler` structure, samplers are executed with the the
commands
```
sample_trajectory!(X, sampler, options=opts);
```
for an in place transformation of `X` or
```
Xvals = sample_trajectory(X₀, sampler);
```
to obtain a full trajectory.  Metropolis methods will also return the running
acceptance rate.

"""
module JuBasicMD

using LinearAlgebra

export sample_trajectory, sample_trajectory!,
    MDOptions,
    RWM, MALA, EM, BBK, ABOBA, BAOAB, HMC, GJF

include("types.jl")
include("sample.jl")
include("utils.jl")
# RWM methods
include("metropolis/zeroth_order/rwm.jl")
# MALA methods
include("metropolis/first_order/mala.jl")
# EM methods
include("nonmetropolis/first_order/em.jl")
# BBK methods
include("nonmetropolis/second_order/bbk.jl")
# ABOBA methods
include("nonmetropolis/second_order/aboba.jl")
# BAOAB methods
include("nonmetropolis/second_order/baoab.jl")
# HMC methods
include("metropolis/second_order/hmc.jl")
# G-JF methods
include("nonmetropolis/second_order/g_jf.jl")
# Verlet methods
# include("nonmetropolis/second_order/verlet.jl")
end
