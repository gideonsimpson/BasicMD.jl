# sample from a harmonic potential in 1D using ABOBA
using Plots
using Printf
using Random
using LinearAlgebra

push!(LOAD_PATH,"../src/")

using JuBasicMD

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

sampler = BAOAB(gradV!, β, γ, M, Δt);

Random.seed!(100);
Q₀ = copy(q₀);
P₀ = copy(p₀);
sample_trajectory!([Q₀, P₀], sampler, options=Options(n_iters=n_iters));
@printf("In Place (Q,P) after %d iterations: (%g,%g)\n",n_iters, Q₀[1], P₀[1]);

Random.seed!(100);
Xvals = sample_trajectory([q₀, p₀], sampler, options=Options(n_iters=n_iters,n_save_iters=n_iters));
Q = Xvals[end][1][1]
P = Xvals[end][2][1]
@printf("(Q,P) after %d iterations: (%g,%g)\n",n_iters, Q[1], P[1]);

Random.seed!(100);
Xvals = sample_trajectory([q₀, p₀], sampler, options=Options(n_iters=n_iters));

histogram([X[1][1] for X in Xvals],label="Samples",normalize=true)
qq=LinRange(-2,2,200)
plot!(qq, sqrt((β)/(2*π))*exp.(-0.5 * β * qq.^2),label="Density")
xlabel!("q")
ylabel!("Frequency")
