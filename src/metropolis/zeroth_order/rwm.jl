struct RWM <: ZerothOrderMetropolisSampler
    V
    β::AbstractFloat
    Δt::AbstractFloat
    gaussian_coeff::AbstractFloat
end

function RWM(V, β, Δt)
    gaussian_coeff = sqrt(2 * Δt /β)
    return RWM(V, β, Δt, gaussian_coeff)
end

mutable struct RWMState{Tx} <:FirstOrderMetropolisSamplerState
    x::Tx
    x_previous::Tx
    V::AbstractFloat
    V_previous::AbstractFloat
    accept::Int
end

function InitState!(initial_x::Tx, sampler::RWM) where Tx

    V = sampler.V(initial_x);
    return RWMState(initial_x, copy(initial_x), V, V, Int(0));
end

function InitState(initial_x::Tx, sampler::RWM) where Tx

    V = sampler.V(initial_x);
    return RWMState(copy(initial_x), copy(initial_x), V, V, Int(0));
end

function UpdateState!(state::RWMState{Tx}, sampler::RWM) where Tx

    @. state.x = state.x_previous +  sampler.gaussian_coeff * randn();
    state.V = sampler.V(state.x);
    a = min(1, exp(sampler.β*(state.V_previous-state.V)))

    if rand()<a
        @. state.x_previous = state.x;
        state.V_previous = state.V;
        state.accept = 1;
    else
        state.accept = 0;
    end
    state
end
