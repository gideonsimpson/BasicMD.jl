# Creating and Using Samplers

```@contents
Pages = ["sample1.md"]
Depth = 4
```

## Sampling Boltzmann 
At their core, all of the included methods sample, approximately, from Botlzmann
distributinos of the type ``\mu(x) \propto e^{-\beta V(x)}``, where the user
must specify the potential, `V(x)`, along with an appropriate inverse temperature, `β`.  

The potential, `V(x)`, must be formulated so that its input argument, `x`
is an array, even if the problem is in ``\mathbb{R}^1``.  For the scalar double well potential, ``V(x) = (x^2-1)^2``, this would be implemneted as:
```
function V(x)
    return (x[1]^2 -1)^2
end
```

### Constructing the Sampler
Having defined the potential, a sampler object must first be defined.  For instance,
```
rwm_sampler = RWM(V, β, Δt);
```
constructs the random walk Metropolis (RWM) sampler for the Boltzmann disribution with ``time step'' `Δt`.

Other samplers require additional arguments.  For instance to use the HMC sampler, we would call
```
sampler = HMC(V, gradV!, β, M, Δt, nΔt);
```
where the additional arguments are:
* `gradV!` -  in-place implementation of ∇V(x), 
* `M` - mass matrix
* `nΔt` - number of Verlet steps of size `Δt` in each HMC iteration.
* 
### Sampling a Trajectory
```@docs 
    sample_trajectory
```


Having created the sampler strucutre and chosen an initial point, `x0`, we call `sample_trajectory`:
```
xvals, avals = sampler_trajectory(x0, sampler);
```
For a Metropolis sampler, like RWM, we return:
* `xvals` - the array of sampled points
* `avals` - the running average of the acceptance rate For non-Metropolis
samplers, the `avals` argument is not returned. 

#### RWM Example

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

#### HMC Example
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

#### Euler-Maruyama Example
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

#### ABOBA Example
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

### Controlling Output
The `sample_trajectory` command will allocates and returns an array of samples.  The number of iterations and the sampling frequency is controlled through the `options` argument.  By default, if this is not specified, the sampler will be called for 10^4 iterations, and record every iteration.  To change this, we construct an [`MDOptions`](@ref) structure and pass that in.

Additionally, one may be in the setting where we do not need to record the full
trajectory, but merely the position at the terminal iterate.  This is handled
with
```@docs 
    sample_trajectory!
```


## Sampling Observables
TBW

## Sampling with Constraints (EXPERIMENTAL)
TBW