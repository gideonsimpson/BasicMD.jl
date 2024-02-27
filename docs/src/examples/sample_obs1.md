# Sampling Observables

```@contents
Pages = ["sample_obs1.md"]
```

These examples are performed using the 2D `EntropicSwitch` potential from [TestLandscapes.jl](https://github.com/gideonsimpson/TestLandscapes.jl)

## RWM Example
```@example
using Plots
using Printf
using Random
using BasicMD
using TestLandscapes

V = x->EntropicSwitch(x);

β = 3.0;
x0 = [-1.0, 0.0];
seed = 100;
Δt = 1e-1;
n_iters = 10^4; # number of samples

sampler = RWM(V, β, Δt);

# construct observables
f₁ = x-> x[1]^2; # second moment of coordinate 1
f₂ = x-> x[2]^2; # second moment of coordinate 2
f₃ = x-> V(x);   # energy
obs = (f₁, f₂, f₃);

Random.seed!(100);

observable_samples = sample_observables(x0, sampler,obs, options=MDOptions(n_iters=n_iters));

plot(1:n_iters, cumsum(observable_samples[1,:])./(1:n_iters), label="E[(X₁)²]")
plot!(1:n_iters, cumsum(observable_samples[2,:])./(1:n_iters), label="E[(X₂)²]")
plot!(1:n_iters, cumsum(observable_samples[3,:])./(1:n_iters), label="E[V(X)]")
xlabel!("Iterate")

```
## HMC Example
```@example
using Plots
using Printf
using Random
using BasicMD
using ForwardDiff
using TestLandscapes

V = x->EntropicSwitch(x);
gradV! = (gradV, x)-> ForwardDiff.gradient!(gradV, V, x);

β = 3.0;
x0 = [-1.0, 0.0];
seed = 100;
M = 1.;
nΔt = 10^1; # number of Verlet steps per HMC iteration
Δt = 1e-1;
n_iters = 10^4; # number of samples

sampler = HMC(V, gradV!, β, M, Δt, nΔt);

# construct observables
f₁ = x-> x[1]^2; # second moment of coordinate 1
f₂ = x-> x[2]^2; # second moment of coordinate 2
f₃ = x-> V(x);   # energy
obs = (f₁, f₂, f₃);

Random.seed!(100);

observable_samples = sample_observables(x0, sampler,obs, options=MDOptions(n_iters=n_iters));

plot(1:n_iters, cumsum(observable_samples[1,:])./(1:n_iters), label="E[(X₁)²]")
plot!(1:n_iters, cumsum(observable_samples[2,:])./(1:n_iters), label="E[(X₂)²]")
plot!(1:n_iters, cumsum(observable_samples[3,:])./(1:n_iters), label="E[V(X)]")
xlabel!("Iterate")

```


## ABOBA Example
!!! note "Coordinates for Inertial Samplers"
    For inertial samplers, which approximate the underdamped Langevin equation, since the trajectory ``x(t) = (q(t), p(t))``, it is essential what when you are evaluating observables, you extract the position or momenum, as appropriate.

```@example
using Plots
using Printf
using Random
using BasicMD
using ForwardDiff

using TestLandscapes

V = x->EntropicSwitch(x);
gradV! = (gradV, x)-> ForwardDiff.gradient!(gradV, V, x);

β = 3.0;
γ = 1.;
M = 1.;
q0 = [-1.0, 0.0];
p0 = [0.0, 0.0]; # as ABOBA is an inertial sampler, we need to specify a momentum
x0 = [copy(q0), copy(p0)];
seed = 100;
Δt = 1e-1;
n_iters = 10^4; # number of samples

sampler = ABOBA(gradV!, β, γ, M, Δt);
opts = MDOptions(n_iters=n_iters);

# construct observables
f₁ = x-> x[1][1]^2; # second moment of coordinate 1
f₂ = x-> x[1][2]^2; # second moment of coordinate 2
f₃ = x-> V(x[1]);   # energy
obs = (f₁, f₂, f₃);

Random.seed!(100);
observable_samples = sample_observables(x0, sampler,obs, options=MDOptions(n_iters=n_iters));

plot(1:n_iters, cumsum(observable_samples[1,:])./(1:n_iters), label="E[(X₁)²]")
plot!(1:n_iters, cumsum(observable_samples[2,:])./(1:n_iters), label="E[(X₂)²]")
plot!(1:n_iters, cumsum(observable_samples[3,:])./(1:n_iters), label="E[V(X)]")
xlabel!("Iterate")
```
