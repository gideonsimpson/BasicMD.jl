
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
* `constraint!` - Modifications of the trajectory i.e., to enforce constraints.
"""
function sample_observables(x₀::Tx, sampler::S, observables::Tuple{Vararg{<:Function,NO}}; options = MDOptions(), constraint! = trivial_constraint!) where {Tx,S<:AbstractSampler,NO}

    state = InitState(x₀, sampler)
    observable_samples = zeros(NO, options.n_save)
    save_index = 1
    for i = 1:options.n_iters
        UpdateState!(state, sampler)
        constraint!(state, i)
        if (mod(i, options.n_save_iters) == 0)
            ntuple(k -> observable_samples[k, save_index] = (observables[k])(state.x), NO)
            # Base.Cartesian.@nexprs $NO k -> observable_samples[k,save_index] = (observables[k])(state.x);
            save_index += 1
        end
    end
    return observable_samples

end

function sample_observables(TO::Type, x₀::Tx, sampler::S, observables::Tuple{Vararg{<:Function,NO}}; options = MDOptions(), constraint! = trivial_constraint!) where {Tx,S<:AbstractSampler,NO}

    state = InitState(x₀, sampler)
    observable_samples = zeros(TO, NO, options.n_save)
    save_index = 1
    for i = 1:options.n_iters
        UpdateState!(state, sampler)
        constraint!(state, i)
        if (mod(i, options.n_save_iters) == 0)
            ntuple(k -> observable_samples[k, save_index] = (observables[k])(state.x), NO)
            # Base.Cartesian.@nexprs $NO k -> observable_samples[k,save_index] = convert(TO,(observables[k])(state.x));
            save_index += 1
        end
    end
    return observable_samples
end

"""
    sample_observables(x₀, sampler, observables; options=MDOptions())

Run the `sampler` starting at `x₀`, evaluating the trajectory on a tuple of
`observables` scalar functions.  Number of iterations and interval between saves
are set using the `options` argument.  Only the computed observables are returned.

### Fields
* `x`         - Starting position for sampler, modified in place
* `sampler`   - Desired sampler
* `recycler`  - Recycling structure for A→B transitions
* `observables` - Observables on which to evaluate the trajectory
### Optional Fields
* `TO`- Observable data type, if needed, should be entered as the first argument
* `options`   - Sampling options, including number of iteration
* `constraint!` - Modifications of the trajectory i.e., to enforce constraints.
"""
function sample_observables(x₀::Tx, sampler::S, recycler::R, observables::Tuple{Vararg{<:Function,NO}}; options = MDOptions(), constraint! = trivial_constraint!) where {Tx,S<:AbstractSampler,R<:AbstractRecycler,NO}

    state = InitState(x₀, sampler)
    observable_samples = zeros(NO, options.n_save)
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
            ntuple(k -> observable_samples[k, save_index] = (observables[k])(state.x), NO)
            # Base.Cartesian.@nexprs $NO k -> observable_samples[k,save_index] = (observables[k])(state.x);
            save_index += 1
        end
    end
    return observable_samples
end

function sample_observables(TO::Type, x₀::Tx, sampler::S, recycler::R, observables::Tuple{Vararg{<:Function,NO}}; options = MDOptions(), constraint! = trivial_constraint!) where {Tx,S<:AbstractSampler,R<:AbstractRecycler,NO}

    state = InitState(x₀, sampler)
    observable_samples = zeros(TO, NO, options.n_save)
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
            ntuple(k -> observable_samples[k, save_index] = (observables[k])(state.x), NO)
            # Base.Cartesian.@nexprs $NO k -> observable_samples[k,save_index] = convert(TO,(observables[k])(state.x));
            save_index += 1
        end
    end
    return observable_samples
end
