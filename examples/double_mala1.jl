# sample from a double well in 1D using HMC
using Plots
using Printf
using Random
using LinearAlgebra
using QuadGK

push!(LOAD_PATH,"../src/")

using JuBasicMD: MALA, MALA!

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
MALA!(X₀, V, gradV!, β, Δt, n_iters);
@printf("In Place X after %d iterations: %g\n",n_iters, X₀[1])

Random.seed!(100);
X, a = MALA(x₀, V,gradV!, β, Δt, n_iters, return_trajectory=false);
@printf("X after %d iterations: %g\n",n_iters, X[1])
@printf("Mean acceptance rate after %d iterations: %g\n",n_iters, a)

Z = quadgk(x->exp(-β*V(x)),-Inf,Inf)[1]
qq = LinRange(-2,2,401)
histogram(Xvals[:],label="Samples",normalize=true)
qq=LinRange(-2,2,200)
plot!(qq, exp.(-β*V.(qq))/Z,label="Density")
xlabel!("x")
ylabel!("Frequency")
