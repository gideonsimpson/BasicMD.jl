# sample from the Muller potential using HMC
using Plots
using Printf
using Random
using LinearAlgebra
using ForwardDiff

push!(LOAD_PATH,"../src/")

using JuBasicMD

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

sampler = HMC(V, gradV!, β, M, Δt, nΔt);

Random.seed!(100);
X₀ = copy(x₀);
sample_trajectory!(X₀, sampler, options=Options(n_iters=n_iters));
@printf("In Place X after %d iterations: (%g, %g)\n",n_iters, X₀[1], X₀[2])

Random.seed!(100);
Xvals, avals = sample_trajectory(x₀, sampler, options=Options(n_iters=n_iters,n_save_iters=n_iters));
X = Xvals[end];
a = avals[end];
@printf("X after %d iterations: (%g, %g)\n",n_iters, X[1], X[2])
@printf("Mean acceptance rate after %d iterations: %g\n",n_iters, a)

Random.seed!(100);
Xvals, avals = sample_trajectory(x₀, sampler, options=Options(n_iters=n_iters));
histogram2d([X[1] for X in Xvals], [X[2] for X in Xvals],normalize=true,color=:viridis)
xlims!(-1.5,1.5)
ylims!(-0.5, 2.0)
xlabel!("x")
ylabel!("y")
