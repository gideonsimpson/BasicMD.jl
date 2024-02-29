struct MALA{TV, TGV, TF<:AbstractFloat} <: FirstOrderMetropolisSampler
    V::TV
    ∇V!::TGV
    β::TF
    Δt::TF
    σ::TF
end

"""
    MALA(V, ∇V!, β, Δt)

Set up the MALA sampler for overdamped Langevin.

### Fields

* V     - Potential
* ∇V!   - In place gradient of the potential
* β     - Inverse temperature
* Δt    - Time step
"""
function MALA(V::TV,∇V!::TGV, β::TF, Δt::TF) where{TV, TGV, TF<:AbstractFloat}
    σ = sqrt(2 * Δt /β)
    return MALA(V, ∇V!, β, Δt, σ)
end

mutable struct MALAState{TF<:AbstractFloat, Tx} <:FirstOrderMetropolisSamplerState
    x::Tx
    x_proposal::Tx
    V::TF
    V_proposal::TF
    ∇V::Tx
    ∇V_proposal::Tx
    accept::Int
end

"""
    MALA_likelihood(X₀, X₁, gradV0, β, Δt)

Compute the likelihood of the MALA proposal X₀→X₁ for
potential with gradient gradV0 at inverse temperature β
and time step Δt
"""
function MALA_likelihood(X₀, X₁, gradV0, β, Δt)
    #norm of increment squared
    norminc2 = 0.0;
    d = length(X₀);
    for i = 1:d
        norminc2 += (X₁[i]-X₀[i] + Δt * gradV0[i])^2;
    end
    return exp( -β * (norminc2) / (4*Δt) );
end

function InitState!(x₀, sampler::MALA)

    V = sampler.V(x₀);
    ∇V = similar(x₀);
    sampler.∇V!(∇V, x₀);
    return MALAState(x₀, copy(x₀),
        V, V, copy(∇V), copy(∇V), Int(0));
end


function InitState(x₀, sampler::MALA)

    V = sampler.V(x₀);
    ∇V = similar(x₀);
    sampler.∇V!(∇V, x₀);
    return MALAState(deepcopy(x₀), similar(x₀),
        V, V, copy(∇V), similar(∇V), Int(0));
end

function UpdateState!(state::MALAState, sampler::MALA)

    @. state.x_proposal = state.x - sampler.Δt * state.∇V + sampler.σ * randn();
    state.V_proposal = sampler.V(state.x_proposal);
    sampler.∇V!(state.∇V_proposal, state.x_proposal);
    g0 = MALA_likelihood(state.x, state.x_proposal, state.∇V, sampler.β, sampler.Δt);
    gp = MALA_likelihood(state.x_proposal, state.x, state.∇V_proposal, sampler.β, sampler.Δt);

    a = min(1, gp/g0 * exp(sampler.β*(state.V-state.V_proposal)));

    if rand()<a
        @. state.x = state.x_proposal;
        @. state.∇V = state.∇V_proposal;
        state.V = state.V_proposal;
        state.accept = 1;
    else
        state.accept = 0;
    end
    state
end
