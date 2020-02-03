# sample from the Muller potential using MALA
using Plots
using Printf
using Random
using LinearAlgebra
using ForwardDiff

push!(LOAD_PATH,"../src/")

using JuBasicMD: MALA, MALA!

β = 0.1;
x₀ = [-0.75, 1.5];
Δt = 1e-3;
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

cfg = ForwardDiff.GradientConfig(V, x₀);
gradV! = (gradV, x)-> ForwardDiff.gradient!(gradV, V, x, cfg);


Random.seed!(100);
X₀ = copy(x₀);
MALA!(X₀, V, gradV!, β, Δt, n_iters);
@printf("In Place X after %d iterations: (%g, %g)\n",n_iters, X₀[1], X₀[2])

Random.seed!(100);
X, a = MALA(x₀, V,gradV!, β, Δt, n_iters, return_trajectory=false);
@printf("X after %d iterations: (%g, %g)\n",n_iters, X[1], X[2])
@printf("Mean acceptance rate after %d iterations: %g\n",n_iters, a)

Random.seed!(100);
Xvals, avals = MALA(x₀, V, gradV!, β, Δt, n_iters, return_trajectory=true);

histogram2d(Xvals[1,:], Xvals[2,:],normalize=true,color=:viridis)
xlims!(-1.5,1.5)
ylims!(-0.5, 2.0)
xlabel!("x")
ylabel!("y")
