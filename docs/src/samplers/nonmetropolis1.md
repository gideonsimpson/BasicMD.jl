# Non-Metropolis Samplers

These methods do __not__ include a Metropolis-Hastings step, and, consequently,
will sample from a distribution, ``\mu_{\Delta t}(x)``, which is a biased
approximation of ``\mu(x) \propto e^{-\beta V(x)}``.   This bias vanishes with
Δt, and is often negligible in comparison to the statistical variance error.

## First Order Methods
These methods are in the spirit of first order in time discretizations.
```@docs
    EM(∇V!::TGV, β::TF, Δt::TF) where{TGV, TF<:AbstractFloat}
```

## Second Order Methods
These methods are in the spirit of second order in time discretizations.
```@docs
    ABOBA(∇V!::TGV, β::TF, γ::TF, M::TM, Δt::TF) where {TGV, TF<:AbstractFloat,TM}
    BAOAB(∇V!::TGV, β::TF, γ::TF, M::TM, Δt::TF) where {TGV, TF<:AbstractFloat,TM}
    BBK(∇V!::TGV, β::TF, γ::TF, M::TM, Δt::TF) where {TGV, TF<:AbstractFloat,TM}
    GJF(∇V!::TGV, β::TF, γ::TF, M::TM, Δt::TF) where {TGV, TF<:AbstractFloat,TM}
```
