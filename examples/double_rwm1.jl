# sample from a double well in 1D using HMC
using Plots
using Printf
using Random
using LinearAlgebra
using ForwardDiff
using QuadGK

push!(LOAD_PATH,"../src/")
using JuBasicMD

include("potentials.jl")

β = 5.0;
x₀ = [-1.0];
seed = 100;
Δt = 1e-1;
n_iters = 10^4; # number of samples

V = x->DoubleWell(x);

sampler = RWM(V, β, Δt);

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
Xvals, avals =sample_trajectory(x₀, sampler, options=MDOptions(n_iters=n_iters));
histogram([X[1] for X in Xvals],label="Samples",normalize=true)
Z = quadgk(x->exp(-β*V(x)),-Inf,Inf)[1]
qq = LinRange(-2,2,401)
plot!(qq, exp.(-β*V.(qq))/Z,label="Density")
xlabel!("x")
ylabel!("Frequency")
