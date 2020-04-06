struct EM{TGV, TF<:AbstractFloat} <: FirstOrderNonMetropolisSampler
    ∇V!::TGV
    β::TF
    Δt::TF
    σ::TF
end

function EM(∇V!::TGV, β::TF, Δt::TF) where{TGV, TF<:AbstractFloat}
    σ = sqrt(2 * Δt /β);
    return EM(∇V!, β, Δt, σ)
end

mutable struct EMState{Tx} <:FirstOrderNonMetropolisSamplerState
    x::Tx
    ∇V::Tx
end

function InitState!(x₀, sampler::EM) where Tx

    ∇V = copy(x₀);
    sampler.∇V!(∇V, x₀);
    return EMState(x₀, copy(∇V));
end

function InitState(x₀, sampler::EM)

    ∇V = copy(x₀);
    sampler.∇V!(∇V, x₀);
    return EMState(copy(x₀), copy(∇V));
end

function UpdateState!(state::EMState, sampler::EM)

    @. state.x = state.x - sampler.Δt * state.∇V + sampler.σ * randn();
    sampler.∇V!(state.∇V, state.x);

    state
end
