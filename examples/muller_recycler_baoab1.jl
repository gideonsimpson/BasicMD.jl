# sample from the Muller potential using BAOAB
using Plots
using Printf
using Random
using LinearAlgebra
using ForwardDiff
using StaticArrays
using BasicMD

include("potentials.jl")

β = 0.1;
γ = 0.1;
M = 1.0;
q₀ = [-0.75, 1.5];
p₀ = [0.0, 0.0];
Δt = 1e-2
n_iters = 10^7;

V = x->Muller(x);
cfg = ForwardDiff.GradientConfig(V, q₀);
gradV! = (gradV, x)-> ForwardDiff.gradient!(gradV, V, x, cfg);

sampler = BAOAB(gradV!, β, γ, M, Δt);

# define the recycling functions
n_recycle = 100;
function inB(state::BasicMD.BAOABState)
    return state.x[1][1]> 0.5;
end
function restartA!(state::BasicMD.BAOABState)
    @. state.x[1] = q₀;
    @. state.x[2] = p₀;
    # gradV!(state.∇V,state.x[1]);
    state
end
recycler = Recycler(restartA!, inB, n_recycle)


Random.seed!(100);
Q₀ = copy(q₀);
P₀ = copy(p₀);
sample_trajectory!([Q₀, P₀], sampler, recycler, options=MDOptions(n_iters=n_iters));
@printf("In Place Q after %d iterations: (%g, %g)\n",n_iters, Q₀[1], Q₀[2])
@printf("In Place P after %d iterations: (%g, %g)\n",n_iters, P₀[1], P₀[2])

Random.seed!(100);
Xvals = sample_trajectory([q₀, p₀], sampler,recycler,  options=MDOptions(n_iters=n_iters,n_save_iters=n_iters));
Q = Xvals[end][1]
P = Xvals[end][2]
@printf("Q after %d iterations: (%g, %g)\n",n_iters, Q[1], Q[2])
@printf("P after %d iterations: (%g, %g)\n",n_iters, P[1], P[2])

Random.seed!(100);
Xvals = sample_trajectory([q₀, p₀], sampler,recycler, options=MDOptions(n_iters=n_iters));
Qvals = [X[1] for X in Xvals];
histogram2d([Q[1] for Q in Qvals], [Q[2] for Q in Qvals],normalize=true,color=:viridis)
xlims!(-1.5,1.5)
ylims!(-0.5, 2.0)
xlabel!("x")
ylabel!("y")
