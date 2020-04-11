struct GJF{TGV, TF<:AbstractFloat, TM} <: SecondOrderNonMetropolisSampler
    ∇V!::TGV
    β::TF
    γ::TF
    M::TM
    Δt::TF
    sqrtM::TM
    a::TF
    b::TF
    σ::TF
end

"""
    GJF(∇V!, β, γ, M, Δt)

Set up the G-JF integrator for inertial Langevin.

### Fields

* ∇V!   - In place gradient of the potential
* β     - Inverse temperature
* γ     - Damping Coefficient
* M     - Mass (either scalar or vector)
* Δt    - Time step
"""
function GJF(∇V!::TGV, β::TF, γ::TF, M::TM, Δt::TF) where {TGV, TF<:AbstractFloat,TM}
    a = (1 - 0.5 * γ * Δt)/(1 + 0.5 * γ * Δt);
    b = 1 / (1 + 0.5 * γ * Δt);
    σ = sqrt(2 * γ * Δt / β);
    sqrtM = sqrt.(M);
    return GJF(∇V!, β, γ, M, Δt, sqrtM, a, b, σ)
end

mutable struct GJFState{Tq, Tx<:AbstractVector{Tq}} <:SecondOrderNonMetropolisSamplerState
    x::Tx
    ∇V::Tq
    ∇V_new::Tq
    ξ::Tq
end

function InitState!(x₀, sampler::GJF)
    ∇V = similar(x₀[1]);
    sampler.∇V!(∇V, x₀[1]);
    return GJFState(x₀, copy(∇V), similar(x₀[1]), similar(x₀[1]));
end

function InitState(x₀, sampler::GJF)

    ∇V = similar(x₀[1]);
    sampler.∇V!(∇V, x₀[1]);
    return GJFState(deepcopy(x₀), copy(∇V), similar(x₀[1]), similar(x₀[1]));
end

function UpdateState!(state::GJFState, sampler::GJF)

    @. state.ξ = randn();

    @. state.x[1] = state.x[1] + sampler.b * sampler.Δt / sampler.M * state.x[2] - 0.5 * sampler.b * sampler.Δt^2 / sampler.M * state.∇V + 0.5 * sampler.b * sampler.Δt / sampler.sqrtM * sampler.σ * state.ξ;

    sampler.∇V!(state.∇V_new, state.x[1]);

    @. state.x[2] = sampler.a * state.x[2] - 0.5 * sampler.Δt * (sampler.a * state.∇V + state.∇V_new) + sampler.b * sampler.sqrtM * sampler.σ * state.ξ;

    @. state.∇V = state.∇V_new;

    state
end
