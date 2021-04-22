# sample from a double well in 1D using ABOBA
using Printf
using Random
using LinearAlgebra
using ForwardDiff
using Statistics
using BasicMD

include("potentials.jl")

β = Float32(5.0);
γ = Float32(1.5);
M = Float32(2.1);
q₀ = Float32.([-1.0]);
p₀ = Float32.([0.0]);
Δt = Float32(1e-1);
n_iters = 10^4;
n_save_iters=10;

V = x->DoubleWell(x);
cfg = ForwardDiff.GradientConfig(V, q₀);
gradV! = (gradV, x)-> ForwardDiff.gradient!(gradV, V, x, cfg);

# second moment and energy obsrevables
f₁ = x-> x[1][1]^2;
f₂ = x-> V(x[1][1]);
observables = (f₁, f₂);

sampler = ABOBA(gradV!, β, γ, M, Δt);

Random.seed!(100);
Xvals = sample_trajectory([q₀, p₀], sampler, options=MDOptions(n_iters=n_iters,n_save_iters=n_save_iters));
f₁_estimate = mean(f₁.(Xvals));
f₂_estimate = mean(f₂.(Xvals));
@printf("f₁ estimate with %d samples: %g\n",length(Xvals), f₁_estimate);
@printf("f₂ estimate with %d samples: %g\n", length(Xvals), f₂_estimate);
Random.seed!(100);
observable_samples = sample_observables(Float32, [q₀, p₀], sampler, observables, options=MDOptions(n_iters=n_iters,n_save_iters=n_save_iters));
@printf("f₁ estimate with %d samples: %g\n",length(observable_samples), mean(observable_samples[1,:]));
@printf("f₂ estimate with %d samples: %g\n", length(observable_samples), mean(observable_samples[2,:]));
