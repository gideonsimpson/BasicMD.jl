module JuBasicMD

using LinearAlgebra
using Printf

"""
    Boltzmann_likelihood(x, V, β)

Compute the unnormalized Boltzmann density, exp(-β V(x))
for potential V at inverse temperature β
"""
function Boltzmann_likelihood(x, V, β)
    w = exp(-β * V(x));
    return w
end

export sample_trajectory, sample_trajectory!,
    Options,
    RWM, MALA, EM, BBK, ABOBA, BAOAB, HMC, Verlet

include("types.jl")
include("sample.jl")
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
# Verlet methods
include("nonmetropolis/second_order/verlet.jl")
end
