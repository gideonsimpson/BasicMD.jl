# sample from a double well in 1D using HMC
using Plots
using Printf
using Random
using LinearAlgebra
using ForwardDiff
using QuadGK
using JuBasicMD

include("potentials.jl")

β = 5.0;
x₀ = [-1.0];
M = 1.0;
seed = 100;
Δt = 2e-1;
n_iters = 10^4; # number of samples
nΔt = 10^1; # number of Verlet steps per HMC iteration

V = x->DoubleWell(x);
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

Z = quadgk(x->exp(-β*V(x)),-Inf,Inf)[1]
qq = LinRange(-2,2,401)
histogram([X[1] for X in Xvals],label="Samples",normalize=true)
plot!(qq, exp.(-β*V.(qq))/Z,label="Density")
xlabel!("x")
ylabel!("Frequency")
