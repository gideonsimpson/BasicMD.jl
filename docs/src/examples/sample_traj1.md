# Sampling Trajectories



```@contents
Pages = ["sample_traj1.md"]
```
Examples of the `sample_trajectory` command with different samplers.

## RWM Example

```@example
using Plots
using Printf
using Random
using BasicMD

function V(x)
    return (x[1]^2 -1)^2
end

β = 5.0;
x0 = [-1.0];
seed = 100;
Δt = 1e-1;
n_iters = 10^4; # number of samples

sampler = RWM(V, β, Δt);

Random.seed!(100);

xvals, avals = sample_trajectory(x0, sampler, options=MDOptions(n_iters=n_iters));

histogram([x_[1] for x_ in xvals],label="RWM Samples",normalize=true)
xlabel!("x")
ylabel!("Frequency")
```

## HMC Example
```@example
using Plots
using Printf
using Random
using BasicMD
using ForwardDiff

function V(x)
    return (x[1]^2 -1)^2
end

gradV! = (gradV, x)-> ForwardDiff.gradient!(gradV, V, x);


β = 5.0;
x0 = [-1.0];
seed = 100;
M = 1.;
nΔt = 10^1; # number of Verlet steps per HMC iteration
Δt = 1e-1;
n_iters = 10^4; # number of samples

sampler = HMC(V, gradV!, β, M, Δt, nΔt);

Random.seed!(100);

xvals, avals = sample_trajectory(x0, sampler, options=MDOptions(n_iters=n_iters));

histogram([x_[1] for x_ in xvals],label="HMC Samples",normalize=true)
xlabel!("x")
ylabel!("Frequency")
```

## Euler-Maruyama Example
```@example
using Plots
using Printf
using Random
using BasicMD
using ForwardDiff

function V(x)
    return (x[1]^2 -1)^2
end

gradV! = (gradV, x)-> ForwardDiff.gradient!(gradV, V, x);


β = 5.0;
x0 = [-1.0];
seed = 100;
Δt = 1e-1;
n_iters = 10^4; # number of samples

sampler = EM(gradV!, β, Δt);

Random.seed!(100);

xvals = sample_trajectory(x0, sampler, options=MDOptions(n_iters=n_iters));

histogram([x_[1] for x_ in xvals],label="EM Samples",normalize=true)
xlabel!("x")
ylabel!("Frequency")
```

## ABOBA Example
```@example
using Plots
using Printf
using Random
using BasicMD
using ForwardDiff

function V(x)
    return (x[1]^2 -1)^2
end

gradV! = (gradV, x)-> ForwardDiff.gradient!(gradV, V, x);


β = 5.0;
γ = 1.;
M = 1.;
q0 = [-1.0];
p0 = [0.0]; # as ABOBA is an inertial sampler, we need to specify a momentum
x0 = [copy(q0), copy(p0)];
seed = 100;
Δt = 1e-1;
n_iters = 10^4; # number of samples

sampler = ABOBA(gradV!, β, γ, M, Δt);
opts = MDOptions(n_iters=n_iters);

Random.seed!(100);
xvals = sample_trajectory(x0, sampler, options= opts);

histogram([x_[1][1] for x_ in xvals],label="ABOBA Samples",normalize=true)
xlabel!("x")
ylabel!("Frequency")
```
