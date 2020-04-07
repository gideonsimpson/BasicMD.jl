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

"""
    HMC(V, ∇V!, β, M, Δt, nΔt)

Set up the HMC sampler for Boltzmann.

### Fields

* V     - Potential
* ∇V!   - In place gradient of the potential
* β     - Inverse temperature
* Δt    - Time step
* nΔt   - Number of time steps to use in each Verlet run
"""
function HMC(V::TV, ∇V!::TGV, β::TF, M::TM, Δt::TF, nΔt::Int) where {TV, TGV, TF<:AbstractFloat,TM}
    sqrtM = sqrt.(M);
    sqrtβ = sqrt(β);
    return HMC(V, ∇V!, β, M, Δt, nΔt, sqrtM,sqrtβ);
end


mutable struct HMCState{Tq, TF<:AbstractFloat} <:SecondOrderMetropolisSamplerState
    x::Tq
    x_proposal::Tq
    V::TF
    V_proposal::TF
    E::TF
    E_proposal::TF
    ∇V::Tq
    p::Tq
    p_mid::Tq
    invMp::Tq
    accept::Int
end

function InitState!(x₀, sampler::HMC)
    V = sampler.V(x₀);
    return HMCState(x₀, copy(x₀), V, V, V, V, similar(x₀), similar(x₀), similar(x₀),similar(x₀), Int(0));
end

function InitState(x₀, sampler::HMC)
    V = sampler.V(x₀);
    return HMCState(copy(x₀), copy(x₀), V, V, V, V, similar(x₀), similar(x₀), similar(x₀),similar(x₀), Int(0));
end

#
# function Verlet!(q, p, ∇V!, M, Δt, nΔt, p_mid, ∇V)
function Verlet!(state::HMCState, sampler::HMC)

    sampler.∇V!(state.∇V, state.x_proposal);
    for _ in 1:sampler.nΔt
        @. state.p_mid = state.p - 0.5 * sampler.Δt * state.∇V;
        @. state.x_proposal = state.x_proposal + sampler.Δt * state.p_mid/sampler.M;
        sampler.∇V!(state.∇V, state.x_proposal);
        @. state.p = state.p_mid - 0.5 * sampler.Δt * state.∇V;
    end
    state
end

function UpdateState!(state::HMCState, sampler::HMC)

    # prepare initial state
    @. state.x_proposal = state.x;
    @. state.p = (sampler.sqrtM/sampler.sqrtβ)  * randn();
    @. state.invMp = state.p / sampler.M;
    state.E = state.V + 0.5 * (state.p ⋅ state.invMp);
    # run Verlet integrator
    Verlet!(state, sampler);
    # compute energy
    state.V_proposal = sampler.V(state.x_proposal);
    @. state.invMp = state.p / sampler.M;
    state.E_proposal = state.V_proposal + 0.5 * (state.p ⋅ state.invMp);
    # accept/reject
    a = min(1, exp( sampler.β *(state.E - state.E_proposal)));
    if rand()<a
        @. state.x = state.x_proposal;
        state.V = state.V_proposal;
        state.accept = 1;
    else
        state.accept = 0;
    end
    state
end
