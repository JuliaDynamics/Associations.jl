export icc, 
    interventional_complexity_causality, 
    InterventionalComplexityCausalityEstimator

"""
    InterventionalComplexityCausalityEstimator

Causality estimators that are based on interventional complexity causality (ICC),
as described in Kathpalia & Nagaraj (2019)[^Kathpalia2019], are subtypes of 
`InterventionalComplexityCausalityEstimator`.

[^Kathpalia2019]: Kathpalia, A., & Nagaraj, N. (2019). Data-based intervention approach for Complexity-Causality measure. PeerJ Computer Science, 5, e196.
"""
abstract type InterventionalComplexityCausalityEstimator end

"""
    icc(x::Vector{Int}, y::Vector{Int}, estimator)

Interventional complexity causality (ICC) is defined as[^Kathpalia2019] *the change in the 
dynamical complexity of time series X when `Xⁱ` is seen to be generated jointly by the 
dynamical evolution of both `Y⁻` and `X⁻` as opposed to by the reality of the dynamical 
evolution of `X⁻` alone"*, where

- `Xⁱ` is a block of current values of `X` (`ΔX` in the original paper), 
- `X⁻` is a block of past values of `X`, not overlapping with and immediately before `Xⁱ`, and 
- `Y⁻` is a block of values from a time series `Y` taken at the same indices as `X⁻`.

Use one of the estimators listed in the documentation.

## Examples

```julia
# Some random binary time series
x = rand(0:1, 500)
y = rand(0:1, 500)
est = CompressionComplexityCausalityEstimator(
    algorithm = EffortToCompress(normalize = true), 
    w = 15, L = 30)
icc(x, y, est)
```

See also: [`CompressionComplexityCausalityEstimator`](@ref), [`EffortToCompress`](@ref)
"""
function interventional_complexity_causality(x, y, estimator::InterventionalComplexityCausalityEstimator) end
const icc = interventional_complexity_causality


