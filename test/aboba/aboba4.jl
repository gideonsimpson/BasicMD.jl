let
    β = 5.0
    q₀ = [-1.0]
    p₀ = [0.0]
    x₀ = [q₀, p₀]
    γ = 1.0
    M = 1.0
    Δt = 1e-1
    
    n_iters = 10^4;
    opts = MDOptions(n_iters=n_iters)

    V = x -> SymmetricDoubleWell(x[1])
    gradV! = (gradV, x) -> ForwardDiff.gradient!(gradV, V, x)

    sampler = ABOBA(gradV!, β, γ, M, Δt)

    f₁ = x -> x[1][1]^2
    f₂ = x -> V(x[1])
    observables = (f₁, f₂)

    Random.seed!(100)
    observable_samples = sample_observables(x₀, sampler, observables, options=opts)
    mean(observable_samples[1, :]) ≈ 0.9365183804236739
end