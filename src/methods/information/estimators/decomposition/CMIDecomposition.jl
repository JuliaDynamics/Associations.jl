export CMIDecomposition

"""
    CMIDecomposition(definition::MultivariateInformationMeasure, 
        est::ConditionalMutualInformationEstimator)

Estimate some multivariate information measure specified by `definition`, by decomposing
it into a combination of conditional mutual information terms. 

## Usage

- Use with [`association`](@ref) to compute a [`MultivariateInformationMeasure`](@ref)
    from input data: [`association`](@ref)`(est::CMIDecomposition, x...)`.
- Use with some [`IndependenceTest`](@ref) to test for independence between variables.


## Description

Each of the conditional mutual information terms are estimated using `est`, which 
can be any [`ConditionalMutualInformationEstimator`](@ref). Finally, these estimates 
are combined according to the relevant decomposition formula.

This estimator is similar to [`EntropyDecomposition`](@ref), but `definition` is expressed as 
conditional mutual information terms instead of entropy terms.

## Examples

- [Example 1](@ref example_TEShannon_CMIDecomposition): Estimating [`TEShannon`](@ref)
    by decomposing it into [`CMIShannon`](@ref) which is estimated using the
    [`FPVP`](@ref) estimator.

See also: [`ConditionalMutualInformationEstimator`](@ref), 
[`MultivariateInformationMeasureEstimator`](@ref),
[`MultivariateInformationMeasure`](@ref).
"""
struct CMIDecomposition{M <: MultivariateInformationMeasure, E} <: DecompositionEstimator{M}
    definition::M # extend API from complexity measures: definition must be the first field of the info estimator.
    est::E # The CMI estimator + measure which `definition` is decomposed into.
end
