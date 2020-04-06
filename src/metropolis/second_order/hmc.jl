struct HMC{TV, TGV, TF<:AbstractFloat, TM} <: SecondOrderMetropolisSampler
    V::TV
    ∇V!::TGV
    β::TF
    M::TM
    Δt::TF
    nΔt::Int
    sqrtM::TM
    sqrtβ::TF
end


function HMC(V::TV, ∇V!::TGV, β::TF, M::TM, Δt::TF, nΔt::Int) where {TV, TGV, TF<:AbstractFloat,TM}
    sqrtM = sqrt.(M);
    sqrtβ = sqrt(β);
    return HMC(V, ∇V!, β, M, Δt, nΔt, sqrtM,sqrtβ);
end


mutable struct HMCState{Tq, TF<:AbstractFloat} <:SecondOrderMetropolisSamplerState
    x::Tq
    x_previous::Tq
    V::TF
    V_previous::TF
    E::TF
    E_previous::TF
    ∇V::Tq
    p::Tq
    p_mid::Tq
    accept::Int
end

function InitState!(x₀, sampler::HMC)
    V = sampler.V(x₀);
    return HMCState(x₀, similar(x₀), V, V, V, V, similar(x₀), similar(x₀), similar(x₀), Int(0));
end

function InitState(x₀, sampler::HMC)
    V = sampler.V(x₀);
    return HMCState(deepcopy(x₀), similar(x₀), V, V, V, V, similar(x₀), similar(x₀), similar(x₀), Int(0));
end

function Verlet!(q::Tq, p::Tq, ∇V!::TGV, M::TM, Δt::TF, nΔt::Int, p_mid::Tq, ∇V::Tq) where {Tq, TGV, TM, TF<:AbstractFloat}

    ∇V!(∇V, q);
    for _ in 1:nΔt
        @. p_mid = p - 0.5 * Δt * ∇V;
        @. q = q + Δt * p_mid/M;
        ∇V!(∇V,q);
        @. p = p_mid - 0.5 * Δt * ∇V;
    end
    q, p
end

function UpdateState!(state::HMCState, sampler::HMC)

    @. state.x = state.x_previous;
    @. state.p = (sampler.sqrtM/sampler.sqrtβ)  * randn();
    state.E_previous = state.V_previous + 0.5 * state.p ⋅ (state.p ./sampler.M);

    # run Verlet integrator
    Verlet!(state.x, state.p, sampler.∇V!, sampler.M, sampler.Δt, sampler.nΔt, state.p_mid, state.∇V);

    # compute energy
    state.V = sampler.V(state.x);
    state.E = state.V + 0.5 * state.p ⋅ (state.p ./sampler.M);

    # accept/reject
    a = min(1, exp( sampler.β *(state.E_previous - state.E)));
    if rand()<a
        @. state.x_previous = state.x;
        state.V_previous = state.V;
        state.accept = 1;
    else
        state.accept = 0;
    end

    state
end
