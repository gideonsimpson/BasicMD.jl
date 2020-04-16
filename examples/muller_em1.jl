# sample from the Muller potential using EM
using Plots
using Printf
using Random
using LinearAlgebra
using ForwardDiff
using StaticArrays

push!(LOAD_PATH,"../src/")

using JuBasicMD

include("potentials.jl")

β = 0.1;
x₀ = [-0.75, 1.5];
Δt = 1e-3;
n_iters = 10^6;

V = x->Muller(x);
cfg = ForwardDiff.GradientConfig(V, x₀);
gradV! = (gradV, x)-> ForwardDiff.gradient!(gradV, V, x, cfg);

sampler = EM(gradV!, β, Δt);

Random.seed!(100);
X₀ = copy(x₀);
sample_trajectory!(X₀, sampler, options=MDOptions(n_iters=n_iters));
@printf("In Place X after %d iterations: (%g, %g)\n",n_iters, X₀[1], X₀[2])

Random.seed!(100);
Xvals = sample_trajectory(x₀, sampler, options=MDOptions(n_iters=n_iters,n_save_iters=n_iters));
X = Xvals[end];
@printf("X after %d iterations: (%g, %g)\n",n_iters, X[1], X[2])

Random.seed!(100);
X_vals = sample_trajectory(x₀, sampler, options=MDOptions(n_iters=n_iters));
histogram2d([X[1] for X in X_vals], [X[2] for X in X_vals],normalize=true,color=:viridis)
xlims!(-1.5,1.5)
ylims!(-0.5, 2.0)
xlabel!("x")
ylabel!("y")
