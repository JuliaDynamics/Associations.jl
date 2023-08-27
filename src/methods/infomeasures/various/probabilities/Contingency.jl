export Contingency

"""
    Contingency <: OutcomeSpace
    Contingency(est::Union{OutcomeSpace, Nothing} = nothing)

`Contingency` is a probabilities estimator that transforms input data to a multidimensional
probability mass function (internally represented as [`ContingencyMatrix`](@ref).

It works directly on raw discrete/categorical data. Alternatively, if a
[`OutcomeSpace`](@ref) `est` for which [`marginal_encodings`](@ref) is implemented
is given, then input data are first discretized before creating the contingency matrix.

!!! note
    The `Contingency` estimator differs from other [`OutcomeSpace`](@ref)s in that
    it's not compatible with [`probabilities`](@ref) and other methods. Instead, it is
    used to construct [`ContingencyMatrix`](@ref), from which probabilities can be
    computed.
"""
Base.@kwdef struct Contingency{E <: Union{Nothing, OutcomeSpace}} <: OutcomeSpace
    est::E = nothing
end
