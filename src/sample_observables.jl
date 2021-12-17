
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
* `constraints` - Constraints on the trajectory
"""
function sample_observables(x₀::Tx, sampler::S, observables::Tuple{Vararg{<:Function,NO}};
    options = MDOptions(), constraints::C = Constraints()) where {Tx,S<:AbstractSampler,NO, C<:AbstractConstraints}

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
            # Base.Cartesian.@nexprs $NO k -> observable_samples[k,save_index] = (observables[k])(state.x);
            save_index += 1
        end
    end
    return observable_samples

end

function sample_observables(TO::Type, x₀::Tx, sampler::S, observables::Tuple{Vararg{<:Function,NO}};
    options = MDOptions(), constraints::C = Constraints()) where {Tx,S<:AbstractSampler,NO, C<:AbstractConstraints}

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
            # Base.Cartesian.@nexprs $NO k -> observable_samples[k,save_index] = convert(TO,(observables[k])(state.x));
            save_index += 1
        end
    end
    return observable_samples
end

