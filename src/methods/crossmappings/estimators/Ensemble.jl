export Ensemble

"""
    Ensemble(; measure::CrossmapMeasure, est::CrossmapEstimator, nreps::Int = 100)

A directive to compute an ensemble analysis, where `measure` (e.g.
[`ConvergentCrossMapping`](@ref)) is computed
using the given estimator `est` (e.g. [`RandomVectors`](@ref))
"""
Base.@kwdef struct Ensemble{M <: CrossmapEstimator{<:CrossmapMeasure}} # todo: perhaps just use a more general Ensemble?
    est::M
    nreps::Int = 100
    function Ensemble(est::M; nreps = 100) where {M}
        new{M}(est, nreps)
    end
end