# Metropolis Samplers

These samplers include a Metropolis-Hastings style step that ensure that the Boltzmann distribution ``\mu(x) \propto e^{-\beta V(x)}`` is _exactly_ targeted; there is no bias associated with, for instance, a finite time step Δt.

## Zeroth Order Methods
These are samplers which do __not__ require the gradienet of the potential, ∇V.

```@docs
    RWM(V::TV, β::TF, Δt::TF) where{TV, TF<:AbstractFloat}
```

## First Order Methods
These are samplers which require the gradienet of the potential, ∇V, and are in
the spirit of first order in time discretizations.

```@docs
    MALA(V::TV,∇V!::TGV, β::TF, Δt::TF) where{TV, TGV, TF<:AbstractFloat}
```

## Second Order Methods
These are samplers which require the gradienet of the potential, ∇V, and are in
the spirit of second order in time discretizations.

```@docs
    HMC(V::TV, ∇V!::TGV, β::TF, M::TM, Δt::TF, nΔt::Int) where {TV, TGV, TF<:AbstractFloat,TM}
```