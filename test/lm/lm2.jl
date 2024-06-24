let
    β = 5.0
    x₀ = [-1.0]
    seed = 100
    Δt = 1e-1
    n_iters = 10^4;
    opts = MDOptions(n_iters=n_iters, n_save_iters=n_iters)

    V = x -> SymmetricDoubleWell(x[1])
    gradV! = (gradV, x) -> ForwardDiff.gradient!(gradV, V, x)

    sampler = LM(gradV!, β, Δt)

    Random.seed!(100)
    Xvals = sample_trajectory(x₀, sampler, options=opts)
    Xvals[end][1] ≈ 1.0408922097822702
end