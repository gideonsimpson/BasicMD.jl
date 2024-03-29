
"""
    sample_observables(x₀, sampler, observables; options=MDOptions())

Run the `sampler` starting at `x₀`, evaluating the trajectory on a tuple of
`observables` scalar functions.  Number of iterations and interval between saves
are set using the `options` argument.  Only the computed observables are returned.

### Fields
* `x`         - Starting position for sampler, modified in place
* `sampler`   - Desired sampler
* `observables` - Observables on which to evaluate the trajectory
### Optional Fields
* `TO`- Observable data type, if needed, should be entered as the first argument
* `options`   - Sampling options, including number of iteration
"""
function sample_observables(x₀::Tx, sampler::S, observables::Tuple{Vararg{<:Function,NO}};
    options = MDOptions()) where {Tx,S<:AbstractSampler,NO}

    state = InitState(x₀, sampler)
    observable_samples = zeros(NO, options.n_save)
    save_index = 1
    for i = 1:options.n_iters
        UpdateState!(state, sampler)
        if (mod(i, options.n_save_iters) == 0)
            ntuple(k -> observable_samples[k, save_index] = (observables[k])(state.x), NO)
            save_index += 1
        end
    end
    return observable_samples

end

function sample_observables(TO::Type, x₀::Tx, sampler::S, observables::Tuple{Vararg{<:Function,NO}};
    options = MDOptions()) where {Tx,S<:AbstractSampler,NO}

    state = InitState(x₀, sampler)
    observable_samples = zeros(TO, NO, options.n_save)
    save_index = 1
    for i = 1:options.n_iters
        UpdateState!(state, sampler)
        if (mod(i, options.n_save_iters) == 0)
            ntuple(k -> observable_samples[k, save_index] = (observables[k])(state.x), NO)
            save_index += 1
        end
    end
    return observable_samples
end

