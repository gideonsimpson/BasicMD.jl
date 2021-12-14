struct MDOptions
    n_iters::Int
    n_save_iters::Int
    n_save::Int
end

"""
    MDOptions(;n_iters = 10^4, n_save_iters=1)

Set options for samplers.

### Fields

* n_iters       - Set the number of iterations of the sampler
* n_save_iters  - Set the frequency at which iterations are saved.  If
                  n_save_iters=1, every iteration is saved.  If n_save_iters=n_iters,
                  only the final iteration is saved.
"""
function MDOptions(; n_iters = 10^4, n_save_iters = 1)

    return MDOptions(n_iters, n_save_iters, floor(Int, n_iters / n_save_iters))
end

"""
    Boltzmann_likelihood(x, V, β)

Compute the unnormalized Boltzmann density, exp(-β V(x))
for potential `V` at inverse temperature `β`

### Fields

* x     - Value at which to evalute the density
* V     - Potential
* β     - Inverse temperature
"""
function Boltzmann_likelihood(x, V, β)
    w = exp(-β * V(x))
    return w
end

"""
`trivial_constraint!` - Trivial constraint function 
### Fields
* `state` - Current state of the sampler
* `i`     - Index of current iterate
"""
function trivial_constraint!(state, i)
    state
end