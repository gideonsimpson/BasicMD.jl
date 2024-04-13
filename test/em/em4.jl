let
    β = 5.0
    x₀ = [-1.0]
    seed = 100
    Δt = 1e-1
    n_iters = 10^4;
    opts = MDOptions(n_iters=n_iters)

    V = x -> SymmetricDoubleWell(x[1])
    gradV! = (gradV, x) -> ForwardDiff.gradient!(gradV, V, x)

    sampler = EM(gradV!, β, Δt)

    f₁ = x -> x[1]^2
    f₂ = x -> V(x)
    observables = (f₁, f₂)

    Random.seed!(100)
    observable_samples = sample_observables(x₀, sampler, observables, options=opts)
    mean(observable_samples[1, :]) ≈ 0.8985651297885703
end