struct Constraints{TA,TB} <: AbstractConstraints
    before_update!::TA
    after_update!::TB
    n_before::Int
    n_after::Int
end

function trivial_constraint!(state::TS) where{TS<:AbstractSamplerState}
    state
end

function Constraints(; before!::TA = trivial_constraint!, after!::TB = trivial_constraint!, n_before::to_indices = 1, n_after::TI = 1) where {TA,TB, TI<:Integer}
    return Constraints(before!, after!, n_before, n_after)
end

