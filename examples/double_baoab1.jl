# sample from a double well in 1D using BAOAB
using Plots
using Printf
using Random
using LinearAlgebra
using QuadGK

push!(LOAD_PATH,"../src/")

using JuBasicMD: BAOAB, BAOAB!

β = 5.0;
γ = 1.5;
M = 2.1;
q₀ = [-1.0];
p₀ = [0.0];
Δt = 1e-1;
n_iters = 10^4;

V = X -> (X⋅X-1)^2;
function gradV!(gradV, X)
    @. gradV = 4 * X * (X^2 -1);
    gradV;
end


Random.seed!(100);
Q₀ = copy(q₀);
P₀ = copy(p₀);
BAOAB!(Q₀, P₀, gradV!, β, γ, M, Δt, n_iters);
@printf("In Place Q after %d iterations: %g\n",n_iters, Q₀[1])

Random.seed!(100);
Q, P = BAOAB(q₀, p₀, gradV!, β, γ, M, Δt, n_iters, return_trajectory=false);
@printf("Q after %d iterations: %g\n",n_iters, Q[1])

Random.seed!(100);
Qvals, Pvals = BAOAB(q₀, p₀, gradV!, β, γ, M, Δt, n_iters);

histogram(Qvals[:],label="Samples",normalize=true)
Z = quadgk(x->exp(-β*V(x)),-Inf,Inf)[1]
qq = LinRange(-2,2,401)
plot!(qq, exp.(-β*V.(qq))/Z,label="Density")
xlabel!("q")
ylabel!("Frequency")
