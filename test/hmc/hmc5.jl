let
    β = 5.0
    x₀ = [-1.0]
    seed = 100
    Δt = 1e-1
    nΔt = 10
    M = 1.0
    n_iters = 10^4
    opts = MDOptions(n_iters=n_iters)

    V = x -> SymmetricDoubleWell(x[1])
    gradV! = (gradV, x) -> ForwardDiff.gradient!(gradV, V, x)


    # define the recycling functions
    a = [-1.0]
    b = 0.5
    function recycler!(state::BasicMD.HMCState)
        if state.x[1] > b
            @. state.x = a
            state.V = V(a)
        end
        state
    end

    recycler = Constraints(recycler!, trivial_constraint!, 1, 1)

    f₁ = x -> x[1]
    f₂ = x -> V(x)
    observables = (f₁, f₂)

    sampler = HMC(V, gradV!, β, M, Δt, nΔt)

    Random.seed!(100)
    observable_samples = sample_observables(x₀, sampler, recycler, observables, options=opts)
    mean(observable_samples[1, :]) ≈ -0.9413081706205281
end