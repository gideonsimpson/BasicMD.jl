# sample from a double well in 1D using HMC
using Plots
using Printf
using Random
using LinearAlgebra
using QuadGK

push!(LOAD_PATH,"../src/")
using JuBasicMD

β = 5.0;
x₀ = [-1.0];
seed = 100;
Δt = 1e-1;
n_iters = 10^4; # number of samples

V = X -> (X⋅X-1)^2;
function gradV!(gradV, X)
    @. gradV = 4 * X * (X^2 -1);
    gradV;
end

Random.seed!(100);
X₀ = copy(x₀);
sampler = MALA(V, gradV!, β, Δt);
sample_trajectory!(X₀, sampler, options=Options(n_iters=n_iters));
@printf("In Place X after %d iterations: %g\n",n_iters, X₀[1])

Random.seed!(100);
X, a = sample_trajectory(x₀, sampler, options=Options(n_iters=n_iters,save_trajectory=false));
@printf("X after %d iterations: %g\n",n_iters, X[1])
@printf("Mean acceptance rate after %d iterations: %g\n",n_iters, a)

Random.seed!(100);
X_vals, a_vals = sample_trajectory(x₀, sampler, options=Options(n_iters=n_iters));
Z = quadgk(x->exp(-β*V(x)),-Inf,Inf)[1]
qq = LinRange(-2,2,401)
histogram([X[1] for X in X_vals],label="Samples",normalize=true)
qq=LinRange(-2,2,200)
plot!(qq, exp.(-β*V.(qq))/Z,label="Density")
xlabel!("x")
ylabel!("Frequency")
