
export MIDecomposition

"""
    MIDecomposition(definition::MultivariateInformationMeasure, 
        est::MutualInformationEstimator)

Estimate some multivariate information measure specified by `definition`, by decomposing
it into a combination of mutual information terms, which are estimated using `est`.

## Usage

- [`information`](@ref)`(est::MIDecomposition, x...)`.

See also: [`MutualInformationEstimator`](@ref), [`MultivariateInformationMeasure`](@ref).
"""
struct MIDecomposition{M <: MultivariateInformationMeasure, E} <: InformationMeasureEstimator{M}
    definition::M # extend API from complexity measures: definition must be the first field of the info estimator.
    est::E # The MI estimator + measure which `definition` is decomposed into.
end