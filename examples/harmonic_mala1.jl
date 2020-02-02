# sample from a harmonic potential in 1D using MALA
using Plots
using Printf
using Random
using LinearAlgebra

push!(LOAD_PATH,"../src/")

using JuBasicMD: MALA, MALA!

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

Random.seed!(100);
X₀ = copy(x₀);
MALA!(X₀, V, gradV!, β, Δt, n_iters);
@printf("In Place X after %d iterations: %g\n",n_iters, X₀[1])

Random.seed!(100);
X, a = MALA(x₀, V,gradV!, β, Δt, n_iters, return_trajectory=false);
@printf("X after %d iterations: %g\n",n_iters, X[1])
@printf("Mean acceptance rate after %d iterations: %g\n",n_iters, a)

Random.seed!(100);
Xvals,avals =MALA(x₀, V, gradV!,β, Δt, n_iters);
histogram(Xvals[:],label="Samples",normalize=true)
qq=LinRange(-2,2,200)
plot!(qq, sqrt((β)/(2*π))*exp.(-0.5 * β * qq.^2),label="Density")
xlabel!("x")
ylabel!("Frequency")
