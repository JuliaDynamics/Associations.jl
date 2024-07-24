export Ensemble

"""
    Ensemble(est::CrossmapEstimator{<:CrossmapMeasure}, nreps::Int = 100)

A directive to compute an ensemble analysis, where `measure` (e.g.
[`ConvergentCrossMapping`](@ref)) is computed
using the given estimator `est` (e.g. [`RandomVectors`](@ref))

## Examples

- [Example 1](@ref example_ConvergentCrossMapping_reproducing_Sugihara_Fig3A): 
    Reproducing Figure 3A from [Sugihara2012](@citet).
"""
Base.@kwdef struct Ensemble{M <: CrossmapEstimator{<:CrossmapMeasure}} # todo: perhaps just use a more general Ensemble?
    est::M
    nreps::Int = 100
    function Ensemble(est::M; nreps = 100) where {M}
        new{M}(est, nreps)
    end
end