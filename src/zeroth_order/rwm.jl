struct RWM <: ZerothOrderSampler
    β::Float64
    Δt::Float64
end

function RWM(β, Δt)
    RWM(β,Δt)
end

mutable struct RWMState{Tx} <:ZerothOrderSamplerState
    x::Tx
    x_previous::Tx
    V_previous::Float64
    accept::Int
end

function InitState(method::RWM, options, V, initial_x::AbstractArray{T}) where T

    initial_V = V(initial_x);
    RWMState(copy(initial_x), copy(initial_x), initial_V, Int(0));
end

function UpdateState!(V, state::RWMState{T})

    @. state.x = state.x_previous sqrt(2 * method.Δt /method.β) * randn();
    state.V = V(state.x)
    V = V(state.x);
    a = min(1, exp(β*(state.V_previous-V)))

    if rand()<a
        @. state.x_previous = state.x

        @. X₀ = Xp;
        V0 = Vp;
    end


end
