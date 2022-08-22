using Statistics
export dynamical_complexity
"""
# Instantaneous dynamical complexity

    dynamical_complexity(Xⁱ, X⁻, algorithm::CompressionComplexityAlgorithm)

Compute the dynamical complexity (DC) of `Xⁱ` given `X⁻` (eq. 1 in Kathpalia & Nagaraj,
2019):

```math
DC(Xⁱ|X⁻) = C(X⁻ \\oplus Xⁱ) - C(X⁻),
```

where ``\\oplus`` signifies vector/array concatenation.

    dynamical_complexity(Xⁱ, X⁻, Y⁻, algorithm::CompressionComplexityCausalityEstimator)

Compute the dynamical complexity of `Xⁱ` given `X⁻` and `Y⁻` (eq. 4 in
Kathpalia & Nagaraj, 2019).

```math
DC(Xⁱ|X⁻, Y⁻) = C(X⁻ \\oplus Xⁱ, Y⁻ \\oplus Xⁱ) - C(X⁻, Y⁻),
```

where ``C(\\cdot, \\cdot)`` is the joint compression complexity.

## Input requirements

Assumes `Xⁱ`, `X⁻` and `Y⁻` are integer-valued sub-segments of pre-symbolized time series
`X` and `Y`, where `X⁻` and `Y⁻` are sampled at the (same) time indices immediately
preceding `Xⁱ`.

## Examples

Estimating dynamical complexity on a single segment using the [`EffortToCompress`](@ref)
algorithm.

```jldoctest; setup = :(using CausalityTools)
using Random; rng = MersenneTwister(1234)
# Some random binary time series segments
x, y = rand(rng, 0:1, 100), rand(rng, 0:1, 100)

# Split segments into past and present
X⁻, Xⁱ = x[1:49], x[50:end]
alg = EffortToCompress(normalize = true)
dynamical_complexity(Xⁱ, X⁻, alg)

# output
-0.09406565656565657
```

# Average dynamical complexity

    dynamical_complexity(f::Function, x, algorithm::CompressionComplexityAlgorithm,
        w::Int, L::Int, step::Int)

Using a constant-width sliding window, compute the instantaneous dynamical complexities of
the integer-valued, pre-symbolized time series `x`, then summarizethe
result using the function `f`. For example, if `f = Statistics.mean`,
the index ``i`` refers to the ``i``-th window, and there are ``N`` windows, we get the
average dynamical complexity

```math
DC(X) = \\dfrac{1}{N}\\sum_{i=1}^{N}DC(Xⁱ|X⁻).
```

The (constant) width of the sliding window is `w + L`, where `L` is the width of the
*present* segment, and `w` is the width of the segment immediately preceding it.

    dynamical_complexity(f::Function, x, y, algorithm::CompressionComplexityAlgorithm,
        w::Int, L::Int, step::Int)

The same as above, but conditioned on a second time series `y`, e.g.

```math
DC(X|Y) = \\dfrac{1}{N}\\sum_{i = 1}^{N}DC(Xⁱ|X⁻, Y⁻),
```

## Examples

Mean dynamical complexity of `x` conditioned on `y`:

```jldoctest; setup = :(using CausalityTools)
using Random; rng = MersenneTwister(1234)
x = rand(rng, 0:1, 200);
y = rand(rng, 0:1, 200);

# Split segments into past and present
alg = EffortToCompress(normalize = true)
w, L, step = 10, 15, 5
dynamical_complexity(Statistics.mean, x, y, alg, w, L, step)

# output
-0.2753968253968254
```
"""
function dynamical_complexity end

function dynamical_complexity(Xⁱ::AbstractVector{J}, X⁻::AbstractVector{J},
        algorithm::CompressionComplexityAlgorithm) where {J <: Integer}
    return compression_complexity([X⁻; Xⁱ], algorithm) -
        compression_complexity(X⁻, algorithm)
end

function dynamical_complexity(Xⁱ::AbstractVector{J}, X⁻::AbstractVector{J}, Y⁻::AbstractVector{J},
    algorithm::CompressionComplexityAlgorithm) where {J <: Integer}
    return compression_complexity([X⁻; Xⁱ], [Y⁻; Xⁱ], algorithm) -
        compression_complexity(X⁻, Y⁻, algorithm)
end

function dynamical_complexity(x, algorithm::CompressionComplexityAlgorithm,
        w::Integer,
        L::Integer,
        step::Integer)
    windows = get_windows(x, w + L, step)
    segment_results = zeros(Float64, length(windows))
    for (i, window) in enumerate(windows)
        current_indices = window[w + 1]:window[w + L]
        past_indices = window[1:w]
        Xⁱ = @views x[current_indices]
        X⁻ = @views x[past_indices]
        segment_results[i] = dynamical_complexity(Xⁱ, X⁻, algorithm)
    end
    return segment_results
end


function dynamical_complexity(x, y,
        algorithm::CompressionComplexityAlgorithm,
        w::Integer,
        L::Integer,
        step::Integer)
    windows = get_windows(x, w + L, step)
    segment_results = zeros(Float64, length(windows))
    for (i, window) in enumerate(windows)
        current_indices = window[w + 1]:window[w + L]
        past_indices = window[1:w]
        Xⁱ = @views x[current_indices]
        X⁻ = @views x[past_indices]
        Y⁻ = @views y[past_indices]
        segment_results[i] = dynamical_complexity(Xⁱ, X⁻, Y⁻, algorithm)
    end
    return segment_results
end

dynamical_complexity(f::Function, x,
    algorithm::CompressionComplexityAlgorithm,
    w::Integer,
    L::Integer,
    step::Integer) =
        f(dynamical_complexity(x, algorithm, w, L, step))

dynamical_complexity(f::Function, x, y,
    algorithm::CompressionComplexityAlgorithm,
    w::Integer,
    L::Integer,
    step::Integer) =
        f(dynamical_complexity(x, y, algorithm, w, L, step))
