
export MIDecomposition

"""
    MIDecomposition(definition::MultivariateInformationMeasure, 
        est::MutualInformationEstimator)

Estimate the [`MultivariateInformationMeasure`](@ref) specified by `definition` by
by decomposing, the measure, if possible, into a combination of mutual information terms.
These terms are individually estimated using the given
[`MutualInformationEstimator`](@ref) `est`, and finally combined to form the final 
value of the measure. 

## Usage

- Use with [`association`](@ref) to compute a [`MultivariateInformationMeasure`](@ref)
    from input data: [`association`](@ref)`(est::MIDecomposition, x...)`.
- Use with some [`IndependenceTest`](@ref) to test for independence between variables.

## Examples

- [Example 1](@ref example_CMIShannon_MIDecomposition): Estimating [`CMIShannon`](@ref)
    using a decomposition into [`MIShannon`](@ref) terms using 
    the [`KraskovStögbauerGrassberger1`](@ref) mutual information estimator.

See also: [`MultivariateInformationMeasureEstimator`](@ref).
"""
struct MIDecomposition{M <: MultivariateInformationMeasure, E} <: DecompositionEstimator{M}
    definition::M # extend API from complexity measures: definition must be the first field of the info estimator.
    est::E # The MI estimator + measure which `definition` is decomposed into.
end
