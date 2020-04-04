
function sample_trajectory!(x::Tx, sampler::S; options=Options()) where {Tx, S<:AbstractSampler}

    state = InitState!(x, sampler);
    for i = 1:options.n_iters
        UpdateState!(state, sampler);
    end
    x
end

function sample_trajectory(initial_x::Tx, sampler::S; options=Options()) where {Tx,  S<:MetropolisSampler}

    n_accept = Int(0);

    state = InitState(initial_x, sampler);

    # allocate memory for samples
    if(options.save_trajectory)
        samples = Tx[similar(initial_x) for j = 1:options.n_iters];
        acceptance_rates = zeros(options.n_iters);
    end

    for i = 1:options.n_iters
        UpdateState!(state,  sampler);
        n_accept+=state.accept;
        if(options.save_trajectory)

            @. samples[i] = state.x;
            acceptance_rates[i] = n_accept/options.n_iters;
        end
    end
    if(options.save_trajectory)
        return samples, acceptance_rates
    else
        return state.x, n_accept/options.n_iters
    end
end
