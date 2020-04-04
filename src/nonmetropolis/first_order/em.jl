struct EM <: FirstOrderNonMetropolisSampler
    ∇V!
    β::AbstractFloat
    Δt::AbstractFloat
    gaussian_coeff::AbstractFloat
end

function EM(∇V!, β, Δt)
    gaussian_coeff = sqrt(2 * Δt /β)
    return EM(∇V!, β, Δt, gaussian_coeff)
end

mutable struct EMState{Tx} <:FirstOrderNonMetropolisSamplerState
    x::Tx
    ∇V::Tx
end

function InitState!(initial_x::Tx, sampler::EM) where Tx

    ∇V = copy(initial_x);
    sampler.∇V!(∇V, initial_x);
    return EMState(initial_x, copy(∇V));
end

function InitState(initial_x::Tx, sampler::EM) where Tx

    ∇V = copy(initial_x);
    sampler.∇V!(∇V, initial_x);
    return EMState(copy(initial_x), copy(∇V));
end

function UpdateState!(state::EMState{Tx}, sampler::EM) where Tx

    @. state.x = state.x - sampler.Δt * state.∇V + sampler.gaussian_coeff * randn();
    sampler.∇V!(state.∇V, state.x);

    state
end
