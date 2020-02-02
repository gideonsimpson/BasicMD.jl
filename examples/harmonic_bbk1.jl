# sample from a harmonic potential in 1D using BBK
using Plots
using Printf
using Random
using LinearAlgebra

push!(LOAD_PATH,"../src/")

using JuBasicMD: BBK, BBK!

β = 5.0;
γ = 1.5;
M = 2.1;
q₀ = [0.0];
p₀ = [0.0];
Δt = 1e-1;
n_iters = 10^4;

function gradV!(gradV, X)
    @. gradV = X;
    gradV;
end

Random.seed!(100);
Q₀ = copy(q₀);
P₀ = copy(p₀);
BBK!(Q₀, P₀, gradV!, β, γ, M, Δt, n_iters);
@printf("In Place Q after %d iterations: %g\n",n_iters, Q₀[1])

Random.seed!(100);
Q, P = BBK(q₀, p₀, gradV!, β, γ, M, Δt, n_iters, return_trajectory=false);
@printf("Q after %d iterations: %g\n",n_iters, Q[1])

Random.seed!(100);
Qvals, Pvals = BBK(q₀, p₀, gradV!, β, γ, M, Δt, n_iters);
histogram(Qvals[:],label="Samples",normalize=true)
qq=LinRange(-2,2,200)
plot!(qq, sqrt((β)/(2*π))*exp.(-0.5 * β * qq.^2),label="Density")
xlabel!("q")
ylabel!("Frequency")
