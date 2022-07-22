# BasicMD.jl Documentation

## Available Samplers
The desired sampler structure must be initialized before any sampling can be performed.

### Metropolis Samplers
```@docs
    RWM(V::TV, β::TF, Δt::TF) where{TV, TF<:AbstractFloat}
    MALA(V::TV,∇V!::TGV, β::TF, Δt::TF) where{TV, TGV, TF<:AbstractFloat}
    HMC(V::TV, ∇V!::TGV, β::TF, M::TM, Δt::TF, nΔt::Int) where {TV, TGV, TF<:AbstractFloat,TM}
```

### Non-Metropolis Samplers
```@docs
    EM(∇V!::TGV, β::TF, Δt::TF) where{TGV, TF<:AbstractFloat}
    ABOBA(∇V!::TGV, β::TF, γ::TF, M::TM, Δt::TF) where {TGV, TF<:AbstractFloat,TM}
    BAOAB(∇V!::TGV, β::TF, γ::TF, M::TM, Δt::TF) where {TGV, TF<:AbstractFloat,TM}
    BBK(∇V!::TGV, β::TF, γ::TF, M::TM, Δt::TF) where {TGV, TF<:AbstractFloat,TM}
    GJF(∇V!::TGV, β::TF, γ::TF, M::TM, Δt::TF) where {TGV, TF<:AbstractFloat,TM}
```

## Options
These options set the number of iterations and the frequency at which data is recorderd.  
```@docs
    MDOptions(; n_iters = 10^4, n_save_iters = 1)
```

## Running Samplers
