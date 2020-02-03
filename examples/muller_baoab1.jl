# sample from the Muller potential using BAOAB
using Plots
using Printf
using Random
using LinearAlgebra
using ForwardDiff

push!(LOAD_PATH,"../src/")

using JuBasicMD: BAOAB, BAOAB!

β = 0.1;
γ = 0.1;
M = [1.5, 0.75];
q₀ = [-0.75, 1.5];
p₀ = [0.0, 0.0];
Δt = 1e-2;
n_iters = 10^6;

function V(x)

    aa = Float64[-1, -1, -6.5, 0.7];
    bb = Float64[0., 0., 11., 0.6];
    cc = Float64[-10., -10., -6.5, 0.7];
    AA = Float64[-200., -100., -170., 15.];
    XX = Float64[1., 0., -0.5, -1.];
    YY = Float64[0., 0.5, 1.5, 1.];

    return ( AA[1]*exp(aa[1]*(x[1]-XX[1])^2+bb[1]*(x[1]-XX[1])*(x[2]-YY[1])+cc[1]*(x[2]-YY[1])^2)
                 +AA[2]*exp(aa[2]*(x[1]-XX[2])^2+bb[2]*(x[1]-XX[2])*(x[2]-YY[2])+cc[2]*(x[2]-YY[2])^2)
                 +AA[3]*exp(aa[3]*(x[1]-XX[3])^2+bb[3]*(x[1]-XX[3])*(x[2]-YY[3])+cc[3]*(x[2]-YY[3])^2)
                 +AA[4]*exp(aa[4]*(x[1]-XX[4])^2+bb[4]*(x[1]-XX[4])*(x[2]-YY[4])+cc[4]*(x[2]-YY[4])^2));
end

cfg = ForwardDiff.GradientConfig(V, q₀);
gradV! = (gradV, x)-> ForwardDiff.gradient!(gradV, V, x, cfg);

Random.seed!(100);
Q₀ = copy(q₀);
P₀ = copy(p₀);
BAOAB!(Q₀, P₀, gradV!, β, γ, M, Δt, n_iters);
@printf("In Place Q after %d iterations: (%g, %g)\n",n_iters, Q₀[1], Q₀[2])

Random.seed!(100);
Q, P = BAOAB(q₀, p₀, gradV!, β, γ, M, Δt, n_iters, return_trajectory=false);
@printf("Q after %d iterations: (%g, %g)\n",n_iters, Q[1], Q[2])

Random.seed!(100);
Qvals, Pvals = BAOAB(q₀, p₀, gradV!, β, γ, M, Δt, n_iters);
histogram2d(Qvals[1,:], Qvals[2,:],normalize=true,color=:viridis)
xlims!(-1.5,1.5)
ylims!(-0.5, 2.0)
xlabel!("x")
ylabel!("y")
