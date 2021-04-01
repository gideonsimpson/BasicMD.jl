# sample from the Muller potential using HMC
using Printf
using Random
using LinearAlgebra
using ForwardDiff
using BasicMD
using Statistics

include("potentials.jl")

β = 0.1;
M = [1.5, 0.75];
x₀ = [-0.75, 1.5];
Δt = 1e-2;
n_iters = 10^6;
nΔt = 10^2;

V = x->Muller(x);
cfg = ForwardDiff.GradientConfig(V, x₀);
gradV! = (gradV, x)-> ForwardDiff.gradient!(gradV, V, x, cfg);

# second moment and energy obsrevables
f₁ = x-> x[1]^2;
f₂ = x-> V(x);
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

