export dynamical_complexity

"""
    dynamical_complexity(Xⁱ::AbstractVector{J}, X⁻::AbstractVector{J}, 
        algorithm::CompressionComplexityAlgorithm)

Compute the dynamical complexity (DC) of `Xⁱ` given `X⁻` (eq. 1 in Kathpalia & Nagaraj,
2019):

```math
DC(Xⁱ|X⁻) = C(X⁻ \\oplus Xⁱ) - C(X⁻),
```

where ``\\oplus`` signifies vector/array concatenation.

    dynamical_complexity(
        Xⁱ::AbstractVector{J}, 
        X⁻::AbstractVector{J}, 
        Y⁻::AbstractVector{J}, 
        algorithm::CompressionComplexityCausalityEstimator)

Compute the dynamical complexity of `Xⁱ` given `X⁻` and `Y⁻` (eq. 4 in 
Kathpalia & Nagaraj, 2019).

```math
DC(Xⁱ|X⁻, Y⁻) = C(X⁻ \\oplus Xⁱ, Y⁻ \\oplus Xⁱ) - C(X⁻, Y⁻),
```

where ``C(\\cdot, \\cdot)`` is the *joint* compression complexity and 
``\\oplus`` signifies vector/array concatenation.

## Input requirements

This function is intended to be used on segments of a longer time series, for example, 
using a sliding window.

Assumes `Xⁱ`, `X⁻` and `Y⁻` are integer-valued sub-segments of pre-symbolized time series
`X` and `Y``, where `X⁻` and `Y⁻` are sampled at the same time indices, from the direct 
past of `Xⁱ`.

## examples

```julia
# Some random binary time series segments
x = rand(0:1, 100)
y = rand(0:1, 100)

# Split segments into past and present
X⁻ = x[1:49]
Xⁱ = x[50:end]
alg = EffortToCompress(normalize = true)
dynamical_complexity(Xⁱ, X⁻, alg)
```
"""
function dynamical_complexity end
const dc = dynamical_complexity

function dynamical_complexity(Xⁱ::AbstractVector{J}, X⁻::AbstractVector{J}, 
        algorithm::CompressionComplexityAlgorithm) where {J <: Integer}
    return compression_complexity([X⁻; Xⁱ], algorithm) - 
        compression_complexity(X⁻, algorithm)
end

function dynamical_complexity(Xⁱ::AbstractVector{J}, X⁻::AbstractVector{J}, Y⁻::AbstractVector{J},
    estimator::CompressionComplexityAlgorithm) where {J <: Integer}
    return compression_complexity([X⁻; Xⁱ], [Y⁻; Xⁱ], algorithm) - 
        compression_complexity(X⁻, Y⁻, algorithm)
end
