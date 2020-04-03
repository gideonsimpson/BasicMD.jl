# sample from a harmonic potential in 1D using RWM
# using Plots
using Printf
using Random
using LinearAlgebra
using Revise
using BenchmarkTools

push!(LOAD_PATH,"../src/")

using JuBasicMD

β = 5.0;
x₀ = [0.0];
seed = 100;
Δt = 1e-1;
n_iters = 10^1;

V = X -> 0.5 * X⋅X;

Random.seed!(100);
X₀ = copy(x₀);

sampler = RWM(V, β, Δt);
sample_trajectory!(X₀, sampler, n_iters);
@printf("In Place X after %d iterations: %g\n",n_iters, X₀[1])

Random.seed!(100);
X_vals, a_vals = sample_trajectory(x₀, sampler, n_iters);

@printf("X after %d iterations: %g\n",n_iters, X_vals[end][1])
@printf("Mean acceptance rate after %d iterations: %g\n",n_iters, a_vals[end])
#
# Random.seed!(100);
# Xvals,avals =RWM(x₀, V, β, Δt, n_iters);
# histogram(Xvals[:],label="Samples",normalize=true)
# qq=LinRange(-2,2,200)
# plot!(qq, sqrt((β)/(2*π))*exp.(-0.5 * β * qq.^2),label="Density")
# xlabel!("x")
# ylabel!("Frequency")
