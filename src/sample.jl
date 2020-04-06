function sample_trajectory!(x::Tx, sampler::S; options=Options()) where {Tx, S<:AbstractSampler}

    state = InitState!(x, sampler);
    for _ in 1:options.n_iters
        UpdateState!(state, sampler);
    end
    x
end

function sample_trajectory(x₀::Tx, sampler::S; options=Options()) where {Tx,  S<:MetropolisSampler}

    n_accept = Int(0);

    state = InitState(x₀, sampler);

    # allocate memory for samples
    samples = Tx[similar(x₀) for i = 1:options.n_save];
    acceptance_rates = zeros(options.n_save);
    save_index = 1;
    for i = 1:options.n_iters
        UpdateState!(state, sampler);
        n_accept+=state.accept;
        if(mod(i,options.n_save_iters)==0)
            @. samples[save_index] = deepcopy(state.x);
            acceptance_rates[save_index] = n_accept/options.n_iters;
            save_index+=1;
        end
    end
    return samples, acceptance_rates
end


function sample_trajectory(x₀::Tx, sampler::S; options=Options()) where {Tx,  S<:NonMetropolisSampler}

    state = InitState(x₀, sampler);
    samples = Tx[similar(x₀) for i = 1:options.n_save];
    save_index = 1;
    for i = 1:options.n_iters
        UpdateState!(state, sampler);
        if(mod(i,options.n_save_iters)==0)
            @. samples[save_index] = deepcopy(state.x);
            save_index+=1;
        end
    end
    return samples
end
