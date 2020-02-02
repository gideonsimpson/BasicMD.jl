# sample from a harmonic potential in 1D using HMC
using Plots
using Printf
using Random
using QuadGK
using LinearAlgebra

push!(LOAD_PATH,"../src/")

using JuBasicMD: HMC

β = 5.0;
x₀ = [0.0];
M = 1.0;
seed = 100;
Δt = 1e-1;
n_iters = 10^4;
nΔt = 10;

V = X -> 0.5 * X⋅X;
function gradV!(gradV, X)
    @. gradV = X;
    gradV;
end


Random.seed!(100);
X, a = HMC(x₀, V, gradV!, β, M, Δt, n_iters,nΔt, return_trajectory=false);
@printf("X after %d iterations: %g\n",n_iters, X[1])
@printf("Mean acceptance rate after %d iterations: %g\n",n_iters, a)

Random.seed!(100);
Xvals,avals =HMC(x₀, V, gradV!, β, M, Δt, n_iters,nΔt, return_trajectory=true);
histogram(Xvals[:],label="Samples",normalize=true)
qq=LinRange(-2,2,200)
plot!(qq, sqrt((β)/(2*π))*exp.(-0.5 * β * qq.^2),label="Density")
xlabel!("x")
ylabel!("Frequency")
