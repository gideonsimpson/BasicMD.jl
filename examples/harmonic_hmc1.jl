# sample from a harmonic potential in 1D using HMC
using Plots
using Printf
using Random
using LinearAlgebra
using ForwardDiff

push!(LOAD_PATH,"../src/")

using JuBasicMD

include("potentials.jl")

β = 5.0;
x₀ = [0.0];
M = 1.0;
seed = 100;
Δt = 1e-1;
n_iters = 10^4;
nΔt = 10;

V = x->Harmonic(x);
cfg = ForwardDiff.GradientConfig(V, x₀);
gradV! = (gradV, x)-> ForwardDiff.gradient!(gradV, V, x, cfg);

sampler = HMC(V, gradV!, β, M, Δt, nΔt);

Random.seed!(100);
X₀ = copy(x₀);
sample_trajectory!(X₀, sampler, options=MDOptions(n_iters=n_iters));
@printf("In Place X after %d iterations: %g\n",n_iters, X₀[1])

Random.seed!(100);
Xvals, avals = sample_trajectory(x₀, sampler, options=MDOptions(n_iters=n_iters,n_save_iters=n_iters));
X = Xvals[end];
a = avals[end];
@printf("X after %d iterations: %g\n",n_iters, X[1])
@printf("Mean acceptance rate after %d iterations: %g\n",n_iters, a)

Random.seed!(100);
Xvals, avals = sample_trajectory(x₀, sampler, options=MDOptions(n_iters=n_iters));
histogram([X[1] for X in Xvals],label="Samples",normalize=true)
qq=LinRange(-2,2,200)
plot!(qq, sqrt((β)/(2*π))*exp.(-0.5 * β * qq.^2),label="Density")
xlabel!("x")
ylabel!("Frequency")
