let 
    β = 5.0
    q₀ = [-1.0]
    p₀ = [0.0]
    x₀ = [q₀, p₀]
    γ = 1.
    M = 1.
    Δt = 1e-1

    n_iters = 10^4
    opts = MDOptions(n_iters=n_iters)

    V = x -> SymmetricDoubleWell(x[1])
    gradV! = (gradV, x) -> ForwardDiff.gradient!(gradV, V, x)

    sampler = BAOAB(gradV!, β, γ, M, Δt)

    Random.seed!(100)
    X₀ = copy(x₀)
    sample_trajectory!(X₀, sampler, options=opts)
    X₀[1][1] ≈ 0.8453671354749337
end