# sample from a double well in 1D using HMC
using Plots
using Printf
using Random
using LinearAlgebra
using QuadGK

push!(LOAD_PATH,"../src/")

using JuBasicMD: RWM, RWM!

β = 5.0;
x₀ = [-1.0];
seed = 100;
Δt = 1e-1;
n_iters = 10^4; # number of samples

function V(X)
    return (X[1]^2 -1)^2
end
#
# function gradV!(gradV, X)
#     @. gradV = 4 * X * (X^2 -1);
#     gradV;
# end

sampler = RWM(V, β, Δt);

Random.seed!(100);
X₀ = copy(x₀);
sample_trajectory!(X₀, sampler, options=Options(n_iters=n_iters));
@printf("In Place X after %d iterations: %g\n",n_iters, X₀[1])

Random.seed!(100);
Xvals, avals = sample_trajectory(x₀, sampler, options=Options(n_iters=n_iters,n_save_iters=n_iters));
X = Xvals[end];
a = avals[end];
@printf("X after %d iterations: %g\n",n_iters, X[1])
@printf("Mean acceptance rate after %d iterations: %g\n",n_iters, a)

Random.seed!(100);
Xvals,avals =RWM(x₀, V, β, Δt, n_iters);

Z = quadgk(x->exp(-β*V(x)),-Inf,Inf)[1]
qq = LinRange(-2,2,401)
histogram(Xvals[:],label="Samples",normalize=true)
qq=LinRange(-2,2,200)
plot!(qq, exp.(-β*V.(qq))/Z,label="Density")
xlabel!("x")
ylabel!("Frequency")
