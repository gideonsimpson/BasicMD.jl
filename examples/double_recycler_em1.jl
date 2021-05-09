# sample from a double well in 1D using HMC
using Plots
using Printf
using Random
using LinearAlgebra
using ForwardDiff
using QuadGK
using BasicMD

include("potentials.jl")

β = 5.0;
x₀ = [-1.0];
seed = 100;
Δt = 1e-2;
n_iters = 10^5; # number of samples

V = x->DoubleWell(x);
cfg = ForwardDiff.GradientConfig(V, x₀);
gradV! = (gradV, x)-> ForwardDiff.gradient!(gradV, V, x, cfg);

sampler = EM(gradV!, β, Δt);

# define the recycling functions
a = [-1.0];
b = 0.9;
function inB(state::BasicMD.EMState)
    return state.x[1]>b
end
function restartA!(state::BasicMD.EMState)
    @. state.x = a;
    state
end


recycler = Recycler(restartA!, inB);

Random.seed!(100);
X₀ = copy(x₀);
sample_trajectory!(X₀, sampler,recycler,  options=MDOptions(n_iters=n_iters));
@printf("In Place X after %d iterations: %g\n",n_iters, X₀[1])

Random.seed!(100);
Xvals = sample_trajectory(x₀, sampler,recycler, options=MDOptions(n_iters=n_iters,n_save_iters=n_iters));
X = Xvals[end];
@printf("X after %d iterations: %g\n",n_iters, X[1])

Random.seed!(100);
X_vals = sample_trajectory(x₀, sampler,recycler, options=MDOptions(n_iters=n_iters));
histogram([X[1] for X in X_vals],label="Samples",normalize=true)
xlabel!("x")
ylabel!("Frequency")

