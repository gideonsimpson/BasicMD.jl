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

V = x -> DoubleWell(x);
cfg = ForwardDiff.GradientConfig(V, x₀);
gradV! = (gradV, x) -> ForwardDiff.gradient!(gradV, V, x, cfg);

sampler = EM(gradV!, β, Δt);

# define the constraint to to keep the trajectory in the set A = [a,b]
a = -1.5;
b = 0.5;
function constraintA!(state::BasicMD.EMState, i)
    if state.x[1] < a
        state.x[1] = a
    end
    if state.x[1] > b
        state.x[1] = b
    end
    state
end

Random.seed!(100);
X₀ = copy(x₀);
sample_trajectory!(X₀, sampler, options = MDOptions(n_iters = n_iters), constraint! = constraintA!);
@printf("In Place X after %d iterations: %g\n", n_iters, X₀[1])

Random.seed!(100);
Xvals = sample_trajectory(x₀, sampler, options = MDOptions(n_iters = n_iters, n_save_iters = n_iters), constraint! = constraintA!);
X = Xvals[end];
@printf("X after %d iterations: %g\n", n_iters, X[1])

Random.seed!(100);
X_vals = sample_trajectory(x₀, sampler, options = MDOptions(n_iters = n_iters),constraint! = constraintA!);
histogram([X[1] for X in X_vals], label = "Samples", normalize = true)
xlabel!("x")
ylabel!("Frequency")

