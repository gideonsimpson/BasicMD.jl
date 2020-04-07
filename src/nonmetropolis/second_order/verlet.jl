struct Verlet{TGV, TF<:AbstractFloat, TM} <: SecondOrderNonMetropolisSampler
    ∇V!::TGV
    M::TM
    Δt::TF
end

function Veret(∇V!::TGV, M::TM, Δt::TF) where {TGV, TF<:AbstractFloat,TM}
    return Verlet(∇V!,M, Δt);
end

mutable struct VerletState{Tq, Tx<:AbstractVector{Tq}, TF<:AbstractFloat} <:SecondOrderNonMetropolisSamplerState
    x::Tx
    ∇V::Tq
    p_mid::Tq
end

function InitState!(x₀, sampler::Verlet)
    ∇V = similar(x₀[1])
    sampler.∇V!(∇V , x₀);
    return VerletState(x₀, copy(∇V) , similar(x₀[1]));
end

function InitState(x₀, sampler::Verlet)
    ∇V = similar(x₀[1])
    sampler.∇V!(∇V , x₀);
    return VerletState(copy(x₀), copy(∇V) , similar(x₀[1]));
end

function UpdateState!(state::VerletState, sampler::Verlet)
    @. state.p_mid = state.p - 0.5 * sampler.Δt * state.∇V;
    # @printf("k = %d: Pmid = %g\n",k, p_mid[1])
    @. state.x = state.x + sampler.Δt * state.p_mid/sampler.M;
    # @printf("k = %d: Qnew = %g\n",k, state.x_proposal[1])
    sampler.∇V!(state.∇V, state.x);
    @. state.p = state.p_mid - 0.5 * sampler.Δt * state.∇V;
    state
end
