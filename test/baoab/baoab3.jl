let
    β = 5.0
    q₀ = [-1.0]
    p₀ = [0.0]
    x₀ = [q₀, p₀]
    γ = 1.0
    M = 1.0
    Δt = 1e-1

    n_iters = 10^4
    opts = MDOptions(n_iters=n_iters)

    V = x -> SymmetricDoubleWell(x[1])
    gradV! = (gradV, x) -> ForwardDiff.gradient!(gradV, V, x)

    # define the recycling functions
    a = [-1.0]
    b = 0.9
    function recycler!(state::BasicMD.BAOABState)
        if state.x[1][1] > b
            @. state.x[1] = a;
            gradV!(state.∇V, a)
        end
        state
    end

    recycler = Constraints(recycler!, trivial_constraint!, 1, 1)

    sampler = BAOAB(gradV!, β, γ, M, Δt)
    Random.seed!(100)
    X₀ = copy(x₀)
    sample_trajectory!(X₀, sampler, recycler, options=opts)

    X₀[1][1] ≈ -1.1695891724498526
end