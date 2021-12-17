struct Constraints{TA,TB} <: AbstractConstraints
    before_update!::TA
    after_update!::TB
    n_before::Int
    n_after::Int
end

"""
    trivial_constraint!

Trival constraint function
"""
function trivial_constraint!(state::TS) where {TS<:AbstractSamplerState}
    state
end

"""
    Constraints

Set constraints on the sampler, if desired
### Fields
* `before! = trivial_constraint!` - Constraint to apply before a step of the sampler
* `after! = trivial_constraint!`  - Constraint to apply after a step of the sampler
* `n_before = 1` - Perform `before!` at `0, n_before-1, 2 n_before-1,...`
* `n_after = 1` - Perform `after!` at `n_after, 2 n_after,...`
"""

function Constraints(; before!::TA = trivial_constraint!, after!::TB = trivial_constraint!, n_before::TI = 1, n_after::TI = 1) where {TA,TB,TI<:Integer}
    return Constraints(before!, after!, n_before, n_after)
end

