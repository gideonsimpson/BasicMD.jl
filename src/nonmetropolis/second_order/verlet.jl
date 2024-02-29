struct Verlet{TGV, TF<:AbstractFloat, TM} <: SecondOrderNonMetropolisSampler
    ∇V!::TGV
    M::TM
    Δt::TF
    
    """
    Verlet(∇V!, M, Δt)
    Set up the Verlet integrator.

    ### Fields

    * ∇V!   - In place gradient of the potential
    * M     - Mass (either scalar or vector)
    * Δt    - Time step
    """
    function Verlet{TGV, TF, TM}(∇V!::TGV, M::TM, Δt::TF) where {TGV, TF<:AbstractFloat, TM} 
        return new(∇V!, M, Δt)
    end 
end

mutable struct VerletState{Tq, Tx<:AbstractVector{Tq}} <:SecondOrderNonMetropolisSamplerState
    x::Tx
    ∇V::Tq
    p_mid::Tq
end

"""
    InitState!(x₀, sampler::Verlet)

TBW
"""
function InitState!(x₀, sampler::Verlet)
    ∇V = similar(x₀[1])
    sampler.∇V!(∇V , x₀[1]);
    return VerletState(x₀, copy(∇V) , similar(x₀[1]));
end

"""
    InitState(x₀, sampler::Verlet)

TBW
"""
function InitState(x₀, sampler::Verlet)
    ∇V = similar(x₀[1])
    sampler.∇V!(∇V , x₀[1]);
    return VerletState(deepcopy(x₀), copy(∇V) , similar(x₀[1]));
end

"""
    UpdateState!(state::VerletState, sampler::Verlet)

TBW
"""
function UpdateState!(state::VerletState, sampler::Verlet)
    @. state.p_mid = state.x[2] - 0.5 * sampler.Δt * state.∇V;
    @. state.x[1] = state.x[1] + sampler.Δt * state.p_mid/sampler.M;
    sampler.∇V!(state.∇V, state.x[1]);
    @. state.x[2] = state.p_mid - 0.5 * sampler.Δt * state.∇V;
    state
end
