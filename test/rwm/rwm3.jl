let
    β = 5.0
    x₀ = [-1.0]
    seed = 100
    Δt = 1e-1
    n_iters = 10^4;
    opts = MDOptions(n_iters=n_iters)

    V = x -> SymmetricDoubleWell(x[1])

    # define the recycling functions
    a = [-1.0]
    b = 0.9
    function recycler!(state::BasicMD.RWMState)
        if state.x[1] > b
            @. state.x = a;
            state.V = V(a);
        end
        state
    end

    recycler = Constraints(recycler!, trivial_constraint!, 1, 1)

    sampler = RWM(V, β, Δt)

    Random.seed!(100)
    X₀ = copy(x₀)
    sample_trajectory!(X₀, sampler, recycler, options=opts)
    X₀[1] ≈ -0.8281049635442811
end