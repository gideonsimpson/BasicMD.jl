struct BBK{TGV, TF<:AbstractFloat, TM} <: SecondOrderNonMetropolisSampler
    ∇V!::TGV
    β::TF
    γ::TF
    M::TM
    Δt::TF
    σ::TF
end

function BBK(∇V!::TGV, β::TF, γ::TF, M::TM, Δt::TF) where {TGV, TF<:AbstractFloat,TM}
    σ = sqrt(γ * Δt / β);
    return BBK(∇V!, β, γ, M, Δt, σ)
end

mutable struct BBKState{Tq, Tx<:AbstractVector{Tq}} <:SecondOrderNonMetropolisSamplerState
    x::Tx
    p_mid::Tq
    ∇V::Tq
end

function InitState!(x₀, sampler::BBK)
    ∇V = similar(x₀[1]);
    sampler.∇V!(∇V, x₀[1]);
    return BBKState(x₀, similar(x₀[1]), copy(∇V));
end

function InitState(x₀, sampler::BBK)

    ∇V = deepcopy(x₀[1]);
    sampler.∇V!(∇V, x₀[1]);
    return BBKState(deepcopy(x₀), similar(x₀[1]), copy(∇V));
end

function UpdateState!(state::BBKState, sampler::BBK)

    @. state.p_mid = state.x[2] - 0.5 * sampler.Δt * state.∇V - 0.5 * sampler.Δt * sampler.γ/sampler.M * state.x[2] + sampler.σ * randn();
    @. state.x[1] = state.x[1] + sampler.Δt/sampler.M * state.p_mid
    sampler.∇V!(state.∇V, state.x[1]);
    @. state.x[2] = (state.p_mid - 0.5 * sampler.Δt * state.∇V + sampler.σ * randn())/(1 + 0.5 * sampler.Δt * sampler.γ/sampler.M);

    state
end
