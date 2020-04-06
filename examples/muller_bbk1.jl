# sample from the Muller potential using BBK
using Plots
using Printf
using Random
using LinearAlgebra
using ForwardDiff
using StaticArrays

push!(LOAD_PATH,"../src/")

using JuBasicMD

β = 0.1;
γ = 0.1;
M = [1.5, 0.75];
q₀ = [-0.75, 1.5];
p₀ = [0.0, 0.0];
Δt = 1e-2;
n_iters = 10^6;

function V(x)

    aa = @SVector [-1, -1, -6.5, 0.7];
    bb = @SVector [0., 0., 11., 0.6];
    cc = @SVector [-10., -10., -6.5, 0.7];
    AA = @SVector [-200., -100., -170., 15.];
    XX = @SVector [1., 0., -0.5, -1.];
    YY = @SVector [0., 0.5, 1.5, 1.];

    return ( AA[1]*exp(aa[1]*(x[1]-XX[1])^2+bb[1]*(x[1]-XX[1])*(x[2]-YY[1])+cc[1]*(x[2]-YY[1])^2)
                 +AA[2]*exp(aa[2]*(x[1]-XX[2])^2+bb[2]*(x[1]-XX[2])*(x[2]-YY[2])+cc[2]*(x[2]-YY[2])^2)
                 +AA[3]*exp(aa[3]*(x[1]-XX[3])^2+bb[3]*(x[1]-XX[3])*(x[2]-YY[3])+cc[3]*(x[2]-YY[3])^2)
                 +AA[4]*exp(aa[4]*(x[1]-XX[4])^2+bb[4]*(x[1]-XX[4])*(x[2]-YY[4])+cc[4]*(x[2]-YY[4])^2));
end

cfg = ForwardDiff.GradientConfig(V, q₀);
gradV! = (gradV, x)-> ForwardDiff.gradient!(gradV, V, x, cfg);

sampler = BBK(gradV!, β, γ, M, Δt)

Random.seed!(100);
Q₀ = copy(q₀);
P₀ = copy(p₀);
sample_trajectory!([Q₀, P₀], sampler, options=Options(n_iters=n_iters));
@printf("In Place Q after %d iterations: (%g, %g)\n",n_iters, Q₀[1], Q₀[2])
@printf("In Place P after %d iterations: (%g, %g)\n",n_iters, P₀[1], P₀[2])

Random.seed!(100);
Xvals = sample_trajectory([q₀, p₀], sampler, options=Options(n_iters=n_iters,n_save_iters=n_iters));
Q = Xvals[end][1]
P = Xvals[end][2]
@printf("Q after %d iterations: (%g, %g)\n",n_iters, Q[1], Q[2])
@printf("P after %d iterations: (%g, %g)\n",n_iters, P[1], P[2])

Random.seed!(100);
Xvals = sample_trajectory([q₀, p₀], sampler, options=Options(n_iters=n_iters));
Qvals = [X[1] for X in Xvals];
histogram2d([Q[1] for Q in Qvals], [Q[2] for Q in Qvals],normalize=true,color=:viridis)
xlims!(-1.5,1.5)
ylims!(-0.5, 2.0)
xlabel!("x")
ylabel!("y")
