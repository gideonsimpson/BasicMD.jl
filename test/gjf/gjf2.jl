let
    β = 5.0
    q₀ = [-1.0]
    p₀ = [0.0]
    x₀ = [q₀, p₀]
    γ = 1.0
    M = 1.0
    Δt = 1e-1

    n_iters = 10^4
    opts = MDOptions(n_iters=n_iters, n_save_iters=n_iters)

    V = x -> SymmetricDoubleWell(x[1])
    gradV! = (gradV, x) -> ForwardDiff.gradient!(gradV, V, x)

    sampler = GJF(gradV!, β, γ, M, Δt)

    Random.seed!(100)
    Xvals = sample_trajectory(x₀, sampler, options=opts)
    Xvals[end][1][1] ≈ 0.8451966064077254
end