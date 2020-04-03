
function sample_trajectory!(x::Tx, sampler::S, n_iters::Int) where {Tx, S<:MetropolisSampler}

    state = InitState(sampler, x);

    for i = 1:n_iters
        UpdateState!(state, sampler);
    end
    @. x = state.x

    x
end

function sample_trajectory(initial_x::Tx, sampler::S, n_iters::Int) where {Tx,  S<:MetropolisSampler}

    n_accept = Int(0);

    state = InitState(sampler, initial_x);

    # allocate memory for samples
    samples = Tx[similar(initial_x) for j = 1:n_iters];
    acceptance_rates = zeros(n_iters);

    for i = 1:n_iters
        UpdateState!(state,  sampler);
        n_accept+=state.accept;
        @. samples[i] = state.x;
        acceptance_rates[i] = n_accept/n_iters;
    end

    return samples, acceptance_rates
end
