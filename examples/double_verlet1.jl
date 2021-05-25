# run a Verlet tarjectory from a double well in 1D
using Plots
using Printf
using LinearAlgebra
using ForwardDiff
using QuadGK
using BasicMD

include("potentials.jl")

M = 1.0;
q₀ = [-0.5];
p₀ = [1.0];
Δt = 1e-1;
n_iters = 10^3;

V = x->DoubleWell(x);
cfg = ForwardDiff.GradientConfig(V, q₀);
gradV! = (gradV, x)-> ForwardDiff.gradient!(gradV, V, x, cfg);

sampler = Verlet(gradV!, M, Δt);

Q₀ = copy(q₀);
P₀ = copy(p₀);
sample_trajectory!([Q₀, P₀], sampler, options=MDOptions(n_iters=n_iters));
@printf("In Place (Q,P) after %d iterations: (%g,%g)\n",n_iters, Q₀[1], P₀[1]);

Xvals = sample_trajectory([q₀, p₀], sampler, 
    options=MDOptions(n_iters=n_iters,n_save_iters=n_iters));
Q = Xvals[end][1][1]
P = Xvals[end][2][1]
@printf("(Q,P) after %d iterations: (%g,%g)\n",n_iters, Q[1], P[1]);

Xvals = sample_trajectory([q₀, p₀], sampler, options=MDOptions(n_iters=n_iters));
Q = Xvals[end][1][1]
P = Xvals[end][2][1]
@printf("(Q,P) after %d iterations: (%g,%g)\n",n_iters, Q[1], P[1]);
plot(Δt *(1:n_iters), [X[1][1] for X in Xvals],label="q")
plot!(Δt *(1:n_iters), [X[2][1] for X in Xvals],label="p")
xlabel!("t")