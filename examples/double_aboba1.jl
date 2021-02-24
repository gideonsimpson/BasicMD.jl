# sample from a double well in 1D using ABOBA
using Plots
using Printf
using Random
using LinearAlgebra
using ForwardDiff
using QuadGK
using JuBasicMD

include("potentials.jl")

β = 5.0;
γ = 1.5;
M = 2.1;
q₀ = [-1.0];
p₀ = [0.0];
Δt = 1e-1;
n_iters = 10^4;

V = x->DoubleWell(x);
cfg = ForwardDiff.GradientConfig(V, q₀);
gradV! = (gradV, x)-> ForwardDiff.gradient!(gradV, V, x, cfg);

sampler = ABOBA(gradV!, β, γ, M, Δt);

Random.seed!(100);
Q₀ = copy(q₀);
P₀ = copy(p₀);
sample_trajectory!([Q₀, P₀], sampler, options=MDOptions(n_iters=n_iters));
@printf("In Place (Q,P) after %d iterations: (%g,%g)\n",n_iters, Q₀[1], P₀[1]);

Random.seed!(100);
Xvals = sample_trajectory([q₀, p₀], sampler, options=MDOptions(n_iters=n_iters,n_save_iters=n_iters));
Q = Xvals[end][1][1]
P = Xvals[end][2][1]
@printf("(Q,P) after %d iterations: (%g,%g)\n",n_iters, Q[1], P[1]);

Random.seed!(100);
Xvals = sample_trajectory([q₀, p₀], sampler, options=MDOptions(n_iters=n_iters));

histogram([X[1][1] for X in Xvals],label="Samples",normalize=true);
Z = quadgk(x->exp(-β*V(x)),-Inf,Inf)[1];
qq = LinRange(-2,2,401)
plot!(qq, exp.(-β*V.(qq))/Z,label="Density")
xlabel!("q")
ylabel!("Frequency")
