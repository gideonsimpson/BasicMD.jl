# sample from a harmonic potential in 1D using RWM
using Plots
using Printf
using Random
using QuadGK
using LinearAlgebra

push!(LOAD_PATH,"../src/")

using JuBasicMD: RWM, RWM!

β = 5.0;
# initial condition is 2D
x₀ = [0.0];
seed = 100;
Δt = 1e-1;
n_iters = 10^4;

V = X -> 0.5 * X⋅X;

Random.seed!(100);
X₀ = copy(x₀);
RWM!(X₀, V, β, Δt, n_iters);
@printf("In Place X after %d iterations: %g\n",n_iters, X₀[1])

Random.seed!(100);
X, a = RWM(x₀, V, β, Δt, n_iters, return_trajectory=false);
@printf("X after %d iterations: %g\n",n_iters, X[1])
@printf("Mean acceptance rate after %d iterations: %g\n",n_iters, a)

Random.seed!(100);
Xvals,avals =RWM(x₀, V, β, Δt, n_iters);
histogram(Xvals[:],label="Samples",normalize=true)
qq=LinRange(-2,2,200)
plot!(qq, sqrt((β)/(2*π))*exp.(-0.5 * β * qq.^2),label="Density")
xlabel!("x")
ylabel!("Frequency")
