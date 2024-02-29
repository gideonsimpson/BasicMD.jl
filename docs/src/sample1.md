# Creating and Using Samplers

```@contents
Pages = ["sample1.md"]
Depth = 4
```

## Sampling Trajectories and Boltzmann
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
    sample_trajectory(::Tx, ::S; options = MDOptions()) where {Tx,S<:BasicMD.MetropolisSampler}
```

Having created the sampler strucutre and chosen an initial point, `x0`, we call `sample_trajectory`:
```
xvals, avals = sampler_trajectory(x0, sampler);
```
For a Metropolis sampler, like RWM, we return:
* `xvals` - the array of sampled points
* `avals` - the running average of the acceptance rate For non-Metropolis
samplers, the `avals` argument is not returned. 

Examples:

```@contents
Pages = ["examples/sample_traj1.md"]
```

### Controlling Output
The `sample_trajectory` command will allocates and returns an array of samples.  The number of iterations and the sampling frequency is controlled through the `options` argument.  By default, if this is not specified, the sampler will be called for 10^4 iterations, and record every iteration.  To change this, we construct an [`MDOptions`](@ref) structure and pass that in.

Additionally, one may be in the setting where we do not need to record the full
trajectory, but merely the position at the terminal iterate.  This is handled
with
```@docs 
    sample_trajectory!(::Tx, ::S; options = MDOptions()) where {Tx,S<:BasicMD.MetropolisSampler}
```


## Sampling Observables

Often, one is not interested in the full trajectory ``\{X_t\}``, but rather the time series observables computed on the trajectory.  Given a function ``f:\mathbb{R}^d\to \mathbb{R}``, it may be satisfactory to simply know ``\{f(X_t)\}``.  When ``d`` is high dimensional, only saving the observavble can cut computational cost and storage. Typically, this is done in order to estimate ergodic averages,
```math
\mathbb{E}_\mu[f(X)] = \int f(x)\mu(dx) \approx \frac{1}{n}\sum_{n=1}^n f(X_{t_n}).
```


This is accomplished using
```@docs 
    sample_observables(x₀::Tx, sampler::S, observables::Tuple{Vararg{<:Function,NO}};
    options = MDOptions()) where {Tx,S<:BasicMD.AbstractSampler,NO}
```

This is quite similar to `sample_trajectory`, except that one must pass an additional argument, a structure containing the desired observables.  Additionally, only the time series of observables is returned, regardless of whether a Metropolis sampler is used or not.  

The `observables` argument should be a tuple of scalar valued functions, i.e. for a single observable:
```
f(x) = x[1]; # first component
obs = (f,);
```
and for multiple observables:
```
f₁(x) = x[1]; # first component
f₂(x) = V(x); # energy
obs = (f₁,f₂);
```
When calling the function,
```
observable_traj = sample_observables(x₀, sampler,obs, options=opts);
```
the returned object, `observable_traj` is a matrix.  Each row corresponding to
an individual observable, recorded at the times specified by the `MDOptions`.

Examples:

```@contents
Pages = ["examples/sample_obs1.md"]
```



## Sampling with Constraints (EXPERIMENTAL)
Elementary tools for enforcing constraints on the system, like ``X_t \in A`` or
``g(X_t)=0`` have been implemented.  This is accomplished by passing a
`constraints` structure to one of `sample_trajectory!`, `sample_trajectory`, or `sample_obsevables`:
```
traj = sample_trajectory(x₀, sampler, constraints);
```
The constraints are constructed with
```@docs
    trivial_constraint!
```

## Recylcing (EXPERIMENTAL)
