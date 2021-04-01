# sample from the Muller potential using BBK
using Printf
using Random
using LinearAlgebra
using ForwardDiff
using BasicMD
using Statistics

include("potentials.jl")

β = 0.1;
γ = 0.1;
M = [1.5, 0.75];
q₀ = [-0.75, 1.5];
p₀ = [0.0, 0.0];
Δt = 1e-2;
n_iters = 10^6;

V = x->Muller(x);
cfg = ForwardDiff.GradientConfig(V, q₀);
gradV! = (gradV, x)-> ForwardDiff.gradient!(gradV, V, x, cfg);

# second moment and energy obsrevables
f₁ = x-> x[1][1]^2;
f₂ = x-> V(x[1]);
observables = (f₁, f₂);

sampler = BBK(gradV!, β, γ, M, Δt)

Random.seed!(100);
Xvals = sample_trajectory([q₀, p₀], sampler, options=MDOptions(n_iters=n_iters,n_save_iters=n_save_iters));
f₁_estimate = mean(f₁.(Xvals));
f₂_estimate = mean(f₂.(Xvals));
@printf("f₁ estimate with %d samples: %g\n",length(Xvals), f₁_estimate);
@printf("f₂ estimate with %d samples: %g\n", length(Xvals), f₂_estimate);
Random.seed!(100);
observable_samples = sample_observables([q₀, p₀], sampler, observables, options=MDOptions(n_iters=n_iters,n_save_iters=n_save_iters));
@printf("f₁ estimate with %d samples: %g\n",length(observable_samples), mean(observable_samples[1,:]));
@printf("f₂ estimate with %d samples: %g\n", length(observable_samples), mean(observable_samples[2,:]));
