# sample from a harmonic potential in 1D using EM
using Plots
using Printf
using Random
using LinearAlgebra

push!(LOAD_PATH,"../src/")

using JuBasicMD

β = 5.0;
x₀ = [0.0];
seed = 100;
Δt = 1e-1;
n_iters = 10^4;

V = X -> 0.5 * X⋅X;
function gradV!(gradV, X)
    @. gradV = X;
    gradV;
end

sampler = EM(gradV!, β, Δt);


Random.seed!(100);
X₀ = copy(x₀);
sample_trajectory!(X₀, sampler, options=Options(n_iters=n_iters));
@printf("In Place X after %d iterations: %g\n",n_iters, X₀[1])

Random.seed!(100);
X = sample_trajectory(x₀, sampler, options=Options(n_iters=n_iters,save_trajectory=false));
@printf("X after %d iterations: %g\n",n_iters, X[1])

Random.seed!(100);
X_vals = sample_trajectory(x₀, sampler, options=Options(n_iters=n_iters));
histogram([X[1] for X in X_vals],label="Samples",normalize=true)
qq=LinRange(-2,2,200)
plot!(qq, sqrt((β)/(2*π))*exp.(-0.5 * β * qq.^2),label="Density")
xlabel!("x")
ylabel!("Frequency")
