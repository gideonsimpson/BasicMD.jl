struct LM{TGV, TF<:AbstractFloat} <: FirstOrderNonMetropolisSampler
    ∇V!::TGV
    β::TF
    Δt::TF
    σ::TF
end

"""
    LM(∇V!, β, Δt)

Set up the LM integrator for overdamped Langevin.

### Fields

* ∇V!   - In place gradient of the potential
* β     - Inverse temperature
* Δt    - Time step
"""
function LM(∇V!::TGV, β::TF, Δt::TF) where{TGV, TF<:AbstractFloat}
    σ = sqrt(2 * Δt /β);
    return LM(∇V!, β, Δt, σ)
end

mutable struct LMState{Tx} <: FirstOrderNonMetropolisSamplerState
    x::Tx
    ξ::Tx
    ξ_new::Tx
    ∇V::Tx
end

function InitState!(x₀, sampler::LM)

    ∇V = similar(x₀);
    ξ = randn(size(x₀));
    ξ_new = similar(x₀);
    sampler.∇V!(∇V, x₀);
    return LMState(x₀, ξ, ξ_new, copy(∇V))
end

function InitState(x₀, sampler::LM)

    ∇V = similar(x₀);
    sampler.∇V!(∇V, x₀);
    ξ = randn(size(x₀));
    ξ_new = similar(x₀)
    return LMState(deepcopy(x₀), ξ,ξ_new, copy(∇V));
end

function UpdateState!(state::LMState, sampler::LM)

    @. state.ξ_new = randn();

    @. state.x = state.x - sampler.Δt * state.∇V + sampler.σ * 0.5 * (state.ξ_new + state.ξ)
    @. state.ξ = state.ξ_new;
    sampler.∇V!(state.∇V, state.x);

    state
end
