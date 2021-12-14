"""
    Recycler(restartA, inB)

Set up the recycler.  This requires that the functions, restartA and inB, take,
as arguments, the same state type as used in the underlying sampling regime

### Fields
* `restartA!`         - Resamples the source state A
* `inB`               - Detects if in state B
* `n_recycle_iters`   - Number of iterations between recycling check
"""
struct Recycler{TA,TB} <: AbstractRecycler
    restartA!::TA
    inB::TB
    n_recycle_iters::Int
end

"""
    Recycler(restartA!::TA, inB::TB; n_recycle_iters = 1)

Construct a recycler.

### Fields
* `restartA!`         - Resamples the source state A
* `inB`               - Detects if in state B
### Optional Fields
* `n_recycle_iters = 1`   - Number of iterations between recycling check
"""
function Recycler(restartA!::TA, inB::TB; n_recycle_iters = 1) where {TA,TB}

    return Recycler(restartA!, inB, n_recycle_iters)
end

"""
    sample_trajectory!(x, sampler, recycler; options=MDOptions())

In place applciation of the `sampler` to `x`.  Number of iterations are
set using the `options` argument.

### Fields

* `x`         - Starting position for sampler, modified in place
* `recycler`  - Recycling structure for A→B transitions
* `sampler`   - Desired sampler
### Optional Fields
* `options`   - Sampling options, including number of iteration
* `constraint!` - Modifications of the trajectory i.e., to enforce constraints
"""
function sample_trajectory!(x::Tx, sampler::S, recycler::R; options = MDOptions(), constraint! = trivial_constraint!) where {Tx,S<:AbstractSampler,R<:AbstractRecycler}

    state = InitState!(x, sampler)
    for i = 1:options.n_iters
        if (mod(i - 1, recycler.n_recycle_iters) == 0)
            if (recycler.inB(state))
                recycler.restartA!(state)
            end
        end
        UpdateState!(state, sampler)
        constraint!(state, i)
    end
    x
end


"""
    sample_trajectory(x₀, sampler, recycler; options=MDOptions())

Run the `sampler` starting at `x₀`.  Number of iterations and interval between
saves are set using the `options` argument.  For Metropolis samplers, the
running acceptance rates are also resturned.


### Fields
* `x`         - Starting position for sampler, modified in place
* `recycler`  - Recycling structure for A→B transitions
* `sampler`   - Desired sampler
* `options`   - Sampling options, including number of iteration
* `constraint!` - Modifications of the trajectory i.e., to enforce constraints
"""
function sample_trajectory(x₀::Tx, sampler::S, recycler::R; options = MDOptions(), constraint! = trivial_constraint!) where {Tx,S<:MetropolisSampler,R<:AbstractRecycler}

    n_accept = Int(0)

    state = InitState(x₀, sampler)

    # allocate memory for samples
    samples = Tx[similar(x₀) for i = 1:options.n_save]
    acceptance_rates = zeros(options.n_save)
    save_index = 1
    for i = 1:options.n_iters
        if (mod(i - 1, recycler.n_recycle_iters) == 0)
            if (recycler.inB(state))
                recycler.restartA!(state)
            end
        end
        UpdateState!(state, sampler)
        constraint!(state, i)
        n_accept += state.accept
        if (mod(i, options.n_save_iters) == 0)
            @. samples[save_index] = deepcopy(state.x)
            acceptance_rates[save_index] = n_accept / i
            save_index += 1
        end
    end
    return samples, acceptance_rates
end


function sample_trajectory(x₀::Tx, sampler::S, recycler::R; options = MDOptions(), constraint! = trivial_constraint!) where {Tx,S<:NonMetropolisSampler,R<:AbstractRecycler}

    state = InitState(x₀, sampler)
    samples = Tx[similar(x₀) for i = 1:options.n_save]
    save_index = 1
    for i = 1:options.n_iters
        if (mod(i - 1, recycler.n_recycle_iters) == 0)
            if (recycler.inB(state))
                recycler.restartA!(state)
            end
        end
        UpdateState!(state, sampler)
        constraint!(state, i)
        if (mod(i, options.n_save_iters) == 0)
            @. samples[save_index] = deepcopy(state.x)
            save_index += 1
        end
    end
    return samples
end
