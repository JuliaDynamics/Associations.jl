export compression_complexity

abstract type CompressionComplexityAlgorithm end

"""
# Regular compression complexity

    compression_complexity(x, algorithm) → N

Compute the compression complexity of the pre-binned/pre-symbolized (integer-valued)
univariate time series `x` using the given `algorithm`.

    compression_complexity(x::Dataset, algorithm, alphabet_size::Int) → N

Multivariate integer time series `x`, given as `Dataset`s, are first transformed into a
1D symbol sequence before computing compression complexity. This transformation assuming
that each variable `xᵢ ∈ x` was symbolized using the same `alphabet_size`. The resulting
1D symbol sequence is (in this implementation) not correct if different alphabet sizes are
used, so ensure during pre-processing that the same alphabet is used (e.g.
`alphabet_size = 2` for a binary time series, and `alphabet_size = 5` for a five-box binned
time series).


# Examples

Quantifying the compression complexity of a univariate symbol sequence using
the [`EffortToCompress`](@ref) algorithm.

```jldoctest; setup = :(using CausalityTools)
compression_complexity([1, 2, 1, 2, 1, 1, 1, 2], EffortToCompress())

# output
5.0
```

Multivariate time series given as [`Dataset`](@ref)s also work, but must be
symbolized using the same alphabet.

```jldoctest; setup = :(using CausalityTools)
x = [1, 2, 1, 2, 1, 1, 1, 2]
y = [2, 1, 2, 2, 2, 1, 1, 2]
alphabet_size = 2
compression_complexity(Dataset(x, y), EffortToCompress(), alphabet_size)

# output
7.0
```

Sliding window estimators automatically handle window creation,
and returns a vector of complexity measures computed on each window.

```jldoctest; setup = :(using CausalityTools)
using Random; rng = MersenneTwister(1234);
x = rand(rng, 0:1, 100)
alg = EffortToCompress(normalize = false)
sw = ConstantWidthSlidingWindow(alg, width = 25, step = 15)
compression_complexity(x, sw)

# output
5-element Vector{Float64}:
 12.0
 14.0
 14.0
 14.0
 11.0
```

# Joint compression complexity

    compression_complexity(x::AbstractVector, y::AbstractVector, algorithm) → N
    compression_complexity(x::Dataset, y::Dataset, algorithm, ax::Int, ay::Int) → N

If a two time series `y` is provided, compute the joint compression
complexity using alphabet size `ax` for `x` and alphabet size `ay` for `y`.

`x` and `y` must be either both integer-valued vectors, or both `Dataset`s (potentially
of different dimensions). If computing the joint compression complexity between a
univariate time series and a multivariate time series, simply wrap the univariate time
series in a `Dataset`.

# Examples

Joint compression complexity between two time series:

```jldoctest; setup = :(using CausalityTools)
using  Random; rng = MersenneTwister(1234);
x = rand(rng, 0:2, 100)
y = rand(rng, 0:2, 100)
alg = EffortToCompress(normalize = true)
compression_complexity(x, y, alg)

# output
0.2222222222222222
```

For multivariate time series, we must specify the alphabet size
for each variable to ensure that results are correct.

```jldoctest; setup = :(using CausalityTools, DelayEmbeddings)
# Multivariate time series X has variables encoded with a 2-letter alphabet,
x1 = [1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0]
x2 = [1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0]
X = Dataset(x1, x2)

# Multivariate time series Y has variables encoded with a 3-letter alphabet
y1 = [1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0]
y2 = [2, 2, 0, 2, 2, 2, 2, 2, 1, 1, 2]
y3 = [2, 2, 0, 2, 2, 2, 1, 2, 1, 1, 2]
Y = Dataset(y1, y2, y3)

alg = EffortToCompress(normalize = true)
compression_complexity(X, Y, alg, 2, 3)

# output
0.9
```

# Returns

See individual algorithms for details on their return values. A `Vector` of compression
complexities is returned if a sliding window algorithm is used - one value per window.
```

See also: [`EffortToCompress`](@ref), [`EffortToCompressSlidingWindow`](@ref).

[^Kathpalia2019]: Kathpalia, A., & Nagaraj, N. (2019). Data-based intervention approach for Complexity-Causality measure. PeerJ Computer Science, 5, e196.
"""
function compression_complexity end
