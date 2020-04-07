struct BAOAB{TGV, TF<:AbstractFloat, TM} <: SecondOrderNonMetropolisSampler
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
    BAOAB(q₀,p₀, gradV!, β, γ, M, Δt, n_iters, return_trajectory=true)

Run the BAOAB integrator for interial Langevin.  If `return_trajectory=true`,
then the entire time series is returned.
"""
function BAOAB(∇V!::TGV, β::TF, γ::TF, M::TM, Δt::TF) where {TGV, TF<:AbstractFloat,TM}

    c₀ = exp(-Δt * γ);
    c₁ = sqrt((1 - exp(-2*γ*Δt))/β);
    sqrtM = sqrt.(M);
    return BAOAB(∇V!, β, γ, M, Δt, c₀, c₁, sqrtM)
end

mutable struct BAOABState{Tq, Tx<:AbstractVector{Tq}} <:SecondOrderNonMetropolisSamplerState
    x::Tx
    q_mid::Tq
    p_mid::Tq
    p̂_mid::Tq
    ∇V::Tq
end

function InitState!(x₀, sampler::BAOAB)
    ∇V = similar(x₀[1]);
    sampler.∇V!(∇V, x₀[1]);
    return BAOABState(x₀, similar(x₀[1]), similar(x₀[1]),similar(x₀[1]), copy(∇V));
end

function InitState(x₀, sampler::BAOAB)

    ∇V = similar(x₀[1]);
    sampler.∇V!(∇V, x₀[1]);
    return BAOABState(deepcopy(x₀), similar(x₀[1]), similar(x₀[1]),similar(x₀[1]), copy(∇V));
end

function UpdateState!(state::BAOABState, sampler::BAOAB)

    @. state.p_mid = state.x[2] - 0.5 * sampler.Δt * state.∇V;
    @. state.q_mid = state.x[1] + 0.5 * sampler.Δt * state.p_mid/sampler.M;
    @. state.p̂_mid = sampler.c₀ * state.p_mid + sampler.c₁ * sampler.sqrtM * randn()
    @. state.x[1] = state.q_mid + 0.5 * sampler.Δt * state.p̂_mid/sampler.M;
    sampler.∇V!(state.∇V,state.x[1]);
    @. state.x[2] = state.p̂_mid - 0.5 * sampler.Δt * state.∇V;

    state
end
