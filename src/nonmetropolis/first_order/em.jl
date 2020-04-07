struct EM{TGV, TF<:AbstractFloat} <: FirstOrderNonMetropolisSampler
    ∇V!::TGV
    β::TF
    Δt::TF
    σ::TF
end

"""
    EM(∇V!, β, γ, M, Δt)

Set up the EM integrator for overdamped Langevin.

### Fields

* ∇V!   - In place gradient of the potential
* β     - Inverse temperature
* Δt    - Time step
"""
function EM(∇V!::TGV, β::TF, Δt::TF) where{TGV, TF<:AbstractFloat}
    σ = sqrt(2 * Δt /β);
    return EM(∇V!, β, Δt, σ)
end

mutable struct EMState{Tx} <:FirstOrderNonMetropolisSamplerState
    x::Tx
    ∇V::Tx
end

function InitState!(x₀, sampler::EM) where Tx

    ∇V = similar(x₀);
    sampler.∇V!(∇V, x₀);
    return EMState(x₀, copy(∇V));
end

function InitState(x₀, sampler::EM)

    ∇V = similar(x₀);
    sampler.∇V!(∇V, x₀);
    return EMState(deepcopy(x₀), copy(∇V));
end

function UpdateState!(state::EMState, sampler::EM)

    @. state.x = state.x - sampler.Δt * state.∇V + sampler.σ * randn();
    sampler.∇V!(state.∇V, state.x);

    state
end
