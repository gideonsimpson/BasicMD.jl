struct RWM{TV, TF<:AbstractFloat} <: ZerothOrderMetropolisSampler
    V::TV
    β::TF
    Δt::TF
    σ::TF
end

function RWM(V::TV, β::TF, Δt::TF) where{TV, TF<:AbstractFloat}
    σ = sqrt(2 * Δt /β);
    return RWM(V, β, Δt, σ)
end

mutable struct RWMState{TF<:AbstractFloat, Tx} <:FirstOrderMetropolisSamplerState
    x::Tx
    x_proposal::Tx
    V::TF
    V_proposal::TF
    accept::Int
end

# function InitState!(x₀::Tx, sampler::RWM{TV, TF}) where {TS, Tx<:AbstractVector{TS},TV, TF<:AbstractFloat}
function InitState!(x₀, sampler::RWM)

    V = sampler.V(x₀);
    return RWMState(x₀, similar(x₀), V, V, Int(0));
end

function InitState(x₀, sampler::RWM)

    V = sampler.V(x₀);
    return RWMState(deepcopy(x₀), similar(x₀), V, V, Int(0));
end

# function UpdateState!(state::RWMState{TF,TS,Tx}, sampler::RWM{TV, TF}) where {TF<:AbstractFloat, TS, Tx<:AbstractVector{TS}, TV}
function UpdateState!(state::RWMState, sampler::RWM)
    @. state.x_proposal = state.x +  sampler.σ * randn();
    state.V_proposal = sampler.V(state.x_proposal);
    a = min(1, exp(sampler.β*(state.V-state.V_proposal)))

    if rand()<a
        @. state.x = state.x_proposal;
        state.V = state.V_proposal;
        state.accept = 1;
    else
        state.accept = 0;
    end
    state
end
