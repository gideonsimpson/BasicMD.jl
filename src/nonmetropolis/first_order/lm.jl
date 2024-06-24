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
    ∇V::Tx
end

function InitState!(x₀, sampler::LM)

    ∇V = similar(x₀);
    ξ = randn(size(x₀));
    sampler.∇V!(∇V, x₀);
    return LMState(x₀, ξ, copy(∇V));
end

function InitState(x₀, sampler::LM)

    ∇V = similar(x₀);
    sampler.∇V!(∇V, x₀);
    ξ = randn(size(x₀))
    return LMState(deepcopy(x₀), ξ, copy(∇V));
end

function UpdateState!(state::LMState, sampler::LM)

    ξ_new = randn(size(state.ξ));

    @. state.x = state.x - sampler.Δt * state.∇V + sampler.σ * .5 *(ξ_new + state.ξ);
    @. state.ξ = ξ_new;
    sampler.∇V!(state.∇V, state.x);

    state
end
