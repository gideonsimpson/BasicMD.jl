let
    β = 5.0
    x₀ = [-1.0]
    seed = 100
    Δt = 1e-1
    nΔt = 10
    M = 1.0
    n_iters = 10^4;
    opts = MDOptions(n_iters=n_iters, n_save_iters=n_iters)

    V = x -> SymmetricDoubleWell(x[1])
    gradV! = (gradV, x) -> ForwardDiff.gradient!(gradV, V, x)

    sampler = HMC(V, gradV!, β, M, Δt, nΔt)

    Random.seed!(100)
    Xvals,_ = sample_trajectory(x₀, sampler, options=opts)
    Xvals[end][1] ≈ 1.0228935878600178
end