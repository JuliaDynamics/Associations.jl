export Contingency

"""
    Contingency <: ProbabilitiesEstimator
    Contingency(est::Union{ProbabilitiesEstimator, Nothing} = nothing)

`Contingency` is a probabilities estimator that transforms input data to a multidimensional
probability mass function (internally represented as [`ContingencyMatrix`](@ref).

It works directly on raw discrete/categorical data. Alternatively, if a
[`ProbabilitiesEstimator`](@ref) `est` for which [`marginal_encodings`](@ref) is implemented
is given, then input data are first discretized before creating the contingency matrix.

!!! note
    `Contingency` estimator differs from other [`ProbabilitiesEstimator`](@ref)s in that
    it's not compatible with [`probabilities`](@ref) and other methods. Instead, it is
    used to construct [`ContingencyMatrix`](@ref), from which probabilities can be
    computed.
"""
Base.@kwdef struct Contingency{E <: Union{Nothing, ProbabilitiesEstimator}} <: ProbabilitiesEstimator
    est::E = nothing
end
