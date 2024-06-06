export CMIDecomposition

"""
    CMIDecomposition(definition::MultivariateInformationMeasure, 
        est::ConditionalMutualInformationEstimator)

Estimate some multivariate information measure specified by `definition`, by decomposing
it into a combination of conditional mutual information terms. Each of these 
terms are then estimated using `est`, which can be any 
[`ConditionalMutualInformationEstimator`](@ref). Finally, these estimates are combined
according to the relevant decomposition formula.

Similar to [`EntropyDecomposition`](@ref), but `definition` is expressed as 
conditional mutual information terms instead of entropy terms.

## Usage

- [`information`](@ref)`(est::CMIDecomposition, x...)`.

## Example 

```julia
using CausalityTools
using Random; rng = MersenneTwister(1234)
x = rand(rng, 100000)
y = rand(rng, 100000) .+ x
z = rand(rng, 100000) .+ y

# Estimate transfer entropy by representing it as a CMI and using the `FPVP` estimator.
est = CMIDecomposition(TEShannon(base = 2), FPVP(k = 3))
information(est, x, z, y) # should be near 0 (and can be negative)
```

See also: [`ConditionalMutualInformationEstimator`](@ref), [`MultivariateInformationMeasure`](@ref).
"""
struct CMIDecomposition{M <: MultivariateInformationMeasure, E} <: DecompositionEstimator{M}
    definition::M # extend API from complexity measures: definition must be the first field of the info estimator.
    est::E # The MI estimator + measure which `definition` is decomposed into.
end
