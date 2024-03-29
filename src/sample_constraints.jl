# sampling functions

"""
    sample_trajectory!(x::Tx, sampler::S, constraints::C; options = MDOptions()) where {Tx,S<:AbstractSampler,C<:AbstractConstraints}

In place applciation of the `sampler` to `x`.  Number of iterations are
set using the `options` argument.

### Fields

* `x`         - Starting position for sampler, modified in place
* `sampler`   - Desired sampler
* `constraints` - Constraints on the trajectory
### Optional Fields
* `options`   - Sampling options, including number of iteration
"""
function sample_trajectory!(x::Tx, sampler::S, constraints::C; options = MDOptions()) where {Tx,S<:AbstractSampler,C<:AbstractConstraints}

    state = InitState!(x, sampler)
    for i = 1:options.n_iters
        if (mod(i - 1, constraints.n_before) == 0)
            constraints.before_update!(state)
        end
        UpdateState!(state, sampler)
        if (mod(i, constraints.n_after) == 0)
            constraints.after_update!(state)
        end
        # constraint!(state, i)
    end
    x
end

"""
    sample_trajectory(x₀::Tx, sampler::S, constraints::C; options = MDOptions()) where {Tx,S<:MetropolisSampler,C<:AbstractConstraints}

Run the `sampler` starting at `x₀`.  Number of iterations and interval between
saves are set using the `options` argument.  For Metropolis samplers, the
running acceptance rates are also resturned.


### Fields
* `x`         - Starting position for sampler
* `sampler`   - Desired sampler
* `constraints` - Constraints on the trajectory
### Optional Fields
* `options`   - Sampling options, including number of iteration
"""
function sample_trajectory(x₀::Tx, sampler::S, constraints::C; options = MDOptions()) where {Tx,S<:MetropolisSampler,C<:AbstractConstraints}

    n_accept = Int(0)

    state = InitState(x₀, sampler)

    # allocate memory for samples
    samples = Tx[similar(x₀) for i = 1:options.n_save]
    acceptance_rates = zeros(options.n_save)
    save_index = 1
    for i = 1:options.n_iters
        if (mod(i - 1, constraints.n_before) == 0)
            constraints.before_update!(state)
        end
        UpdateState!(state, sampler)
        if (mod(i, constraints.n_after) == 0)
            constraints.after_update!(state)
        end
        n_accept += state.accept
        if (mod(i, options.n_save_iters) == 0)
            @. samples[save_index] = deepcopy(state.x)
            acceptance_rates[save_index] = n_accept / i
            save_index += 1
        end
    end
    return samples, acceptance_rates
end

function sample_trajectory(x₀::Tx, sampler::S, constraints::C; options = MDOptions()) where {Tx,S<:NonMetropolisSampler,C<:AbstractConstraints}

    state = InitState(x₀, sampler)
    samples = Tx[similar(x₀) for i = 1:options.n_save]
    save_index = 1
    for i = 1:options.n_iters
        if (mod(i - 1, constraints.n_before) == 0)
            constraints.before_update!(state)
        end
        UpdateState!(state, sampler)
        if (mod(i, constraints.n_after) == 0)
            constraints.after_update!(state)
        end
        if (mod(i, options.n_save_iters) == 0)
            @. samples[save_index] = deepcopy(state.x)
            save_index += 1
        end
    end
    return samples
end

"""
    sample_observables(x₀::Tx, sampler::S, constraints::C, observables::Tuple{Vararg{<:Function,NO}};
    options = MDOptions()) where {Tx,S<:AbstractSampler,NO,C<:AbstractConstraints}

Run the `sampler` starting at `x₀`, evaluating the trajectory on a tuple of
`observables` scalar functions.  Number of iterations and interval between saves
are set using the `options` argument.  Only the computed observables are returned.

### Fields
* `x`         - Starting position for sampler, modified in place
* `sampler`   - Desired sampler
* `constraints` - Constraints on the trajectory
* `observables` - Observables on which to evaluate the trajectory
### Optional Fields
* `TO`- Observable data type, if needed, should be entered as the first argument
* `options`   - Sampling options, including number of iteration
"""
function sample_observables(x₀::Tx, sampler::S, constraints::C, observables::Tuple{Vararg{<:Function,NO}};
    options = MDOptions()) where {Tx,S<:AbstractSampler,NO,C<:AbstractConstraints}

    state = InitState(x₀, sampler)
    observable_samples = zeros(NO, options.n_save)
    save_index = 1
    for i = 1:options.n_iters
        if (mod(i - 1, constraints.n_before) == 0)
            constraints.before_update!(state)
        end
        UpdateState!(state, sampler)
        if (mod(i, constraints.n_after) == 0)
            constraints.after_update!(state)
        end
        if (mod(i, options.n_save_iters) == 0)
            ntuple(k -> observable_samples[k, save_index] = (observables[k])(state.x), NO)
            save_index += 1
        end
    end
    return observable_samples

end

function sample_observables(TO::Type, x₀::Tx, sampler::S, constraints::C, observables::Tuple{Vararg{<:Function,NO}};
    options = MDOptions()) where {Tx,S<:AbstractSampler,NO,C<:AbstractConstraints}

    state = InitState(x₀, sampler)
    observable_samples = zeros(TO, NO, options.n_save)
    save_index = 1
    for i = 1:options.n_iters
        if (mod(i - 1, constraints.n_before) == 0)
            constraints.before_update!(state)
        end
        UpdateState!(state, sampler)
        if (mod(i, constraints.n_after) == 0)
            constraints.after_update!(state)
        end
        if (mod(i, options.n_save_iters) == 0)
            ntuple(k -> observable_samples[k, save_index] = (observables[k])(state.x), NO)
            save_index += 1
        end
    end
    return observable_samples
end

