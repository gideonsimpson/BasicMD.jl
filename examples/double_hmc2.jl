# sample from a double well in 1D using HMC
using Printf
using Random
using LinearAlgebra
using ForwardDiff
using QuadGK
using BasicMD
using Statistics

include("potentials.jl")

β = 5.0;
x₀ = [-1.0];
M = 1.0;
seed = 100;
Δt = 2e-1;
n_iters = 10^4; # number of samples
nΔt = 10^1; # number of Verlet steps per HMC iteration
n_save_iters=10;

V = x->DoubleWell(x);
cfg = ForwardDiff.GradientConfig(V, x₀);
gradV! = (gradV, x)-> ForwardDiff.gradient!(gradV, V, x, cfg);

# second moment and energy obsrevables
f₁ = x-> x[1]^2;
f₂ = x-> V(x[1]);
observables = (f₁, f₂);

sampler = HMC(V, gradV!, β, M, Δt, nΔt);

Random.seed!(100);
Xvals,_ = sample_trajectory(x₀, sampler, options=MDOptions(n_iters=n_iters,n_save_iters=n_save_iters));
f₁_estimate = mean(f₁.(Xvals));
f₂_estimate = mean(f₂.(Xvals));
@printf("f₁ estimate with %d samples: %g\n",length(Xvals), f₁_estimate);
@printf("f₂ estimate with %d samples: %g\n", length(Xvals), f₂_estimate);
Random.seed!(100);
observable_samples = sample_observables(x₀, sampler, observables, options=MDOptions(n_iters=n_iters,n_save_iters=n_save_iters));
@printf("f₁ estimate with %d samples: %g\n",length(observable_samples), mean(observable_samples[1,:]));
@printf("f₂ estimate with %d samples: %g\n", length(observable_samples), mean(observable_samples[2,:]));
