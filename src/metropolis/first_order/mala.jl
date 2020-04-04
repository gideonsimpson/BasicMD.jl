struct  MALA <: FirstOrderMetropolisSampler
    V
    ∇V!
    β::AbstractFloat
    Δt::AbstractFloat
    gaussian_coeff::AbstractFloat
end


function MALA(V, ∇V!, β, Δt)
    gaussian_coeff = sqrt(2 * Δt /β)
    return MALA(V, ∇V!, β, Δt, gaussian_coeff)
end

mutable struct MALAState{Tx} <:ZerothOrderMetropolisSamplerState
    x::Tx
    x_previous::Tx
    V::AbstractFloat
    V_previous::AbstractFloat
    ∇V::Tx
    ∇V_previous::Tx
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
    for i = 1:length(X₀)
        norminc2 += (X₁[i]-X₀[i] + Δt * gradV0[i])^2;
    end
    return exp( -β * (norminc2) / (4*Δt) );
end

function InitState!(initial_x::Tx, sampler::MALA) where Tx

    V = sampler.V(initial_x);
    ∇V = copy(initial_x);
    sampler.∇V!(∇V, initial_x);
    return MALAState(initial_x, copy(initial_x),
        V, V, copy(∇V), copy(∇V), Int(0));
end


function InitState(initial_x::Tx, sampler::MALA) where Tx

    V = sampler.V(initial_x);
    ∇V = copy(initial_x);
    sampler.∇V!(∇V, initial_x);
    return MALAState(copy(initial_x), copy(initial_x),
        V, V, copy(∇V), copy(∇V), Int(0));
end

function UpdateState!(state::MALAState{Tx}, sampler::MALA) where Tx

    @. state.x = state.x_previous - sampler.Δt * state.∇V_previous + sampler.gaussian_coeff * randn();
    state.V = sampler.V(state.x);
    sampler.∇V!(state.∇V,state.x);
    g0 = MALA_likelihood(state.x_previous, state.x, state.∇V_previous, sampler.β, sampler.Δt);
    gp = MALA_likelihood(state.x, state.x_previous, state.∇V, sampler.β, sampler.Δt);

    a = min(1, gp/g0 * exp(sampler.β*(state.V_previous-state.V)));

    if rand()<a
        @. state.x_previous = state.x;
        @. state.∇V_previous = state.∇V;
        state.V_previous = state.V;
        state.accept = 1;
    else
        state.accept = 0;
    end
    state
end
