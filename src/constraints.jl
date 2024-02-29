"""
Constraints{TA, TB, TI}(before_update!::TA, after_update!::TB, n_before::TI, n_after::TI) where {TA, TB, TI<:Integer}

Set constraints on the sampler, if desired
### Fields
* `before_update!` - Constraint to apply before a step of the sampler
* `after_update!`  - Constraint to apply after a step of the sampler
* `n_before` - Perform `before!` at `0, n_before-1, 2 n_before-1,...`
* `n_after` - Perform `after!` at `n_after, 2 n_after,...`
"""
struct Constraints{TA,TB,TI<:Integer} <: AbstractConstraints
    before_update!::TA
    after_update!::TB
    n_before::TI
    n_after::TI
end

"""
    trivial_constraint!(state::TS) where {TS<:AbstractSamplerState}

Trival constraint function
"""
function trivial_constraint!(state::TS) where {TS<:AbstractSamplerState}
    state
end
