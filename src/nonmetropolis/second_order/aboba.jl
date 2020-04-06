struct ABOBA{TGV, TF<:AbstractFloat, TM} <: SecondOrderNonMetropolisSampler
    ∇V!::TGV
    β::TF
    γ::TF
    M::TM
    Δt::TF
    c₀::TF
    c₁::TF
    sqrtM::TM
end

"""
    ABOBA(q₀,p₀, gradV!, β, γ, M, Δt, n_iters, return_trajectory=true)

Run the ABOBA integrator for interial Langevin.  If `return_trajectory=true`,
then the entire time series is returned.
"""
function ABOBA(∇V!::TGV, β::TF, γ::TF, M::TM, Δt::TF) where {TGV, TF<:AbstractFloat,TM}

    c₀ = exp(-Δt * γ);
    c₁ = sqrt((1 - exp(-2*γ*Δt))/β);
    sqrtM = sqrt.(M);
    return ABOBA(∇V!, β, γ, M, Δt, c₀, c₁, sqrtM)
end

mutable struct ABOBAState{Tq, Tx<:AbstractVector{Tq}} <:SecondOrderNonMetropolisSamplerState
    x::Tx
    q_mid::Tq
    p_mid::Tq
    p̂_mid::Tq
    ∇V_mid::Tq
end

function InitState!(x₀, sampler::ABOBA)
    return ABOBAState(x₀, similar(x₀[1]), similar(x₀[1]), similar(x₀[1]), similar(x₀[1]));
end

function InitState(x₀, sampler::ABOBA)
    return ABOBAState(deepcopy(x₀), similar(x₀[1]), similar(x₀[1]),similar(x₀[1]), similar(x₀[1]));
end

function UpdateState!(state::ABOBAState, sampler::ABOBA)

    @. state.q_mid = state.x[1] + 0.5 * sampler.Δt * state.x[2]/sampler.M;
    sampler.∇V!(state.∇V_mid,state.q_mid);
    @. state.p_mid = state.x[2] - 0.5 * sampler.Δt * state.∇V_mid;
    @. state.p̂_mid = sampler.c₀ * state.p_mid + sampler.c₁ * sampler.sqrtM * randn()
    @. state.x[2] = state.p̂_mid - 0.5 * sampler.Δt * state.∇V_mid;
    @. state.x[1] = state.q_mid + 0.5 * sampler.Δt * state.x[2]/sampler.M;

    state
end
