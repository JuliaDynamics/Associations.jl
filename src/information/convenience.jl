export condmutualinfo

"""
    condmutualinfo(est::CondiitionalMutualInformationEstimator, x, y, z) → cmi::Real
    condmutualinfo(est::MIDecomposition, x, y, z) → cmi::Real
    condmutualinfo(est::JointProbabilities, x, y, z) → cmi::Real
    condmutualinfo(est::EntropyDecomposition, x, y, z) → cmi::Real

Estimate some [`ConditionalMutualInformation`](@ref) `*`, ``CMI_*(X, Y | Z)``, using the 
given `estimator`
type ``*`` using  the estimator `est`, where type of the CMI is controlled by `est.definition`.

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

| Estimator                      | Principle                    | [`CMIShannon`](@ref) | [`CMIRenyiPoczos`](@ref) | [`CMIRenyiJizba`](@ref) | [`CMIRenyiSarbu`](@ref) | [`CMITsallisPapapetrou`](@ref) |
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

# Discrete estimation from joint probabilities
information(JointProbabilities(CMIShannon(), ValueBinning(3)), x, z, y)

# Four-entropies decomposition with `Kraskov` differential entropy estimator for all terms.
information(EntropyDecomposition(CMIShannon(), Kraskov(k = 20)), x, z, y)

# Two-mutual-informations decomposition with `GaoKannanOhViswanath` MI estimator for both terms.
information(MIDecomposition(CMIShannon(), GaoKannanOhViswanath(k = 20)), x, z, y)
```
"""
function condmutualinfo(est::MultivariateInformationMeasureEstimator{<:ConditionalMutualInformation}, x, y, z)
    return information(est, x, y, z)
end