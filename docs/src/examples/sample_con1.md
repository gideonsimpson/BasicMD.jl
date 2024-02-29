# Sampling with Constraints

```@contents
Pages = ["sample_con1.md"]
```
## Doublewell Nonequilibrium Steady State
For the double well potential, ``V(x) = (x^2-1)^2``, consider the overdamped
Langevin dynamics:
```math
dX_t = - \nabla V(X_t)dt + \sqrt{2\beta^{-1}}dW_t.
```
with ``X_0 = -1``.  Now, we wish to find the nonequiblrium steady state (NESS)
such that when the process arrives in the set ``B=[0.5, \infty)``, it restarts
at ``-1``.  We can accomplish this with an Euler-Maruyama discretization and
imposing the constraint to reset the problem resulting in the modified process
``\tilde{X}_t``.  The following code plots the distribution.  This can also be
used to find the MFPT.
```@example
using Plots
using Printf
using Random
using ForwardDiff
using BasicMD

function V(x)
    return (x[1]^2 -1)^2
end

gradV! = (gradV, x)-> ForwardDiff.gradient!(gradV, V, x);

β = 5.0;
x₀ = [-1.0];
seed = 100;
Δt = 1e-2;
n_iters = 10^4; # number of samples

sampler = EM(gradV!, β, Δt);

# define the recycling function and the constraint
a = [-1.0];
b = 0.5;
function recycler!(state::BasicMD.EMState)
    if state.x[1] > b
        @. state.x = a
        gradV!(state.∇V, a)
    end
    state
end

recycler = Constraints(recycler!, trivial_constraint!, 1, 1);

Random.seed!(100);
X_vals = sample_trajectory(x₀, sampler, recycler, options = MDOptions(n_iters = n_iters));
histogram([X[1] for X in X_vals], label = "Samples", normalize = true, bins = 25)
xlabel!("x")
ylabel!("Frequency")
title!("NESS")

```

