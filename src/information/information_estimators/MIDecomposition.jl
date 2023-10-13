
export MIDecomposition

"""
    MIDecomposition(definition::MultivariateInformationMeasure, 
        est::MutualInformationEstimator)

Estimate some multivariate information measure specified by `definition`, by decomposing
it into a combination of mutual information terms, which are estimated using `est`,
which can be any [`MutualInformationEstimator`](@ref).

Similar to [`EntropyDecomposition`](@ref), but `definition` is expressed as 
mutual information terms instead of entropy terms.

## Usage

- [`information`](@ref)`(est::MIDecomposition, x...)`.

## Example 

```julia
using CausalityTools
using Random; rng = MersenneTwister(1234)
x = rand(rng, 100000)
y = rand(rng, 100000) .+ x
z = rand(rng, 100000) .+ y

est = MIDecomposition(CMIShannon(), KSG1(MIShannon(base = 2), k = 3))
information(est, x, z, y) # should be near 0 (and can be negative)
```

See also: [`MutualInformationEstimator`](@ref), [`MultivariateInformationMeasure`](@ref).
"""
struct MIDecomposition{M <: MultivariateInformationMeasure, E} <: MultivariateInformationMeasureEstimator{M}
    definition::M # extend API from complexity measures: definition must be the first field of the info estimator.
    est::E # The MI estimator + measure which `definition` is decomposed into.
end