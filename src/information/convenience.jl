export conditional_entropy
export mutualinfo
export condmutualinfo

function conditional_entropy(est::MultivariateInformationMeasureEstimator{<:ConditionalEntropy}, x, y)
    return information(est, x, y)
end

function mutualinfo(est::MultivariateInformationMeasureEstimator{<:MutualInformation}, x, y)
    return information(est, x, y)
end

"""
    condmutualinfo(est::ConditionalMutualInformationEstimator, x, y, z) → cmi::Real

Estimate the conditional mutual information (CMI) measure specified by `est.definition`
between `x` and `y`, given `z`.

## Estimation

The following estimators are dedicated [`ConditionalMutualInformationEstimator`](@ref).

| Estimator                    | Principle         | [`CMIShannon`](@ref) | [`CMIRenyiPoczos`](@ref) |
| ---------------------------- | ----------------- | :------------------: | :----------------------: |
| [`FPVP`](@ref)               | Nearest neighbors |          ✓          |            x             |
| [`MesnerShalizi`](@ref)      | Nearest neighbors |          ✓          |            x             |
| [`Rahimzamani`](@ref)        | Nearest neighbors |          ✓          |            x             |
| [`PoczosSchneiderCMI`](@ref) | Nearest neighbors |          x           |            ✓            |
| [`GaussianCMI`](@ref)        | Parametric        |          ✓          |            x             |

The CMI measures may also be estimated using any of the following generic estimators:

| Estimator                      | Principle                    | [`CMIShannon`](@ref) | [`CMIRenyiPoczos`](@ref) | [`CMIRenyiJizba`](@ref) | [`CMIRenyiSarbu`](@ref) | [`CMITsallis`](@ref) |
| ------------------------------ | ---------------------------- | :------------------: | :----------------------: | :---------------------: | :---------------------: | :------------------: |
| [`JointProbabilities`](@ref)   | Discrete joint pmf           |          ✓          |            x             |           ✓            |           ✓            |          ✓          |
| [`EntropyDecomposition`](@ref) | Four-entropies decomposition |          ✓          |            x             |           ✓            |            x            |          x           |
| [`MIDecomposition`](@ref)      | Two-MI decomposition         |          ✓          |            x             |            x            |            x            |          x           |

## Examples

```julia
using CausalityTools
using Random; rng = MersenneTwister(1234)
x = rand(rng, 1000)
y = rand(rng, 1000) .+ x
z = rand(rng, 1000) .+ y

# Dedicated estimators
condmutualinfo(FPVP(CMIShannon(base = 2), k = 20), x, z, y)
condmutualinfo(PoczosSchneiderCMI(CMIRenyiPoczos(), k = 20), x, z, y)

# Generic estimators
information(JointProbabilities(CMIShannon(), ValueBinning(3)), x, z, y)

```
"""
function condmutualinfo(est::MultivariateInformationMeasureEstimator{<:ConditionalMutualInformation}, x, y, z)
    return information(est, x, y, z)
end