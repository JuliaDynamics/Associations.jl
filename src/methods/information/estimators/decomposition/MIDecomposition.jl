
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
    from input data.

## Examples

One common application is computing Shannon-type conditional mutual information.
It can be decomposed as a sum of mutual information terms, which we can each 
estimate with any dedicated [`MutualInformationEstimator`](@ref) estimator.

```julia
using CausalityTools
using Random; rng = MersenneTwister(1234)
x = rand(rng, 100000)
y = rand(rng, 100000) .+ x
z = rand(rng, 100000) .+ y

est = MIDecomposition(CMIShannon(), KSG1(MIShannon(base = 2), k = 3))
association(est, x, z, y) # should be near 0 (and can be negative)
```

See also: [`EntropyDecomposition`](@ref).
"""
struct MIDecomposition{M <: MultivariateInformationMeasure, E} <: DecompositionEstimator{M}
    definition::M # extend API from complexity measures: definition must be the first field of the info estimator.
    est::E # The MI estimator + measure which `definition` is decomposed into.
end
