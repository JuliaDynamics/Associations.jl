
"""
    estimate(e::MultivariateInformationMeasure, est::InformationMeasureEstimator,
        x::VectorOrStateSpaceSet...)

Given some input data `x`, estimate some multivariate information measure using the given
[`InformationMeasureEstimator`](@ref) `est`, with respect to the generalized entropy `e`.
"""
function estimate(measure::MultivariateInformationMeasure, args...; kwargs...) end
