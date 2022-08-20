export compression_complexity

abstract type CompressionComplexityAlgorithm end

"""
# Regular compression complexity

    compression_complexity(x, algorithm) → N

Compute the compression complexity of the pre-binned/pre-symbolized (integer-valued)
univariate (`Vector{Int}`) or multivariate (`Dataset`) time series `x`, using the given 
`algorithm`. For multivariate time series, the algorithm's `alphabet_size` must be a 
nonzero integer giving the number of possible states each state vector can take.

# Joint compression complexity

    compression_complexity(x, y, algorithm) → N
    compression_complexity(x::Dataset, y::Dataset, algorithm, ax::Int, ay::Int) → N

If a second time series `y` is provided, compute the joint compression
complexity using alphabet size `ax` for `x` and alphabet size `ay` for `y`.

# Returns

See individual algorithms for details on their return values.
    
A `Vector` of compression complexities is returned if a sliding window algorithm is used - 
one value per window.

# Examples

Let's use the effort-to-compress [`EffortToCompress`](@ref) algorithm to quantify the 
compression complexity of a univariate symbol sequence.

```julia
# A univariate time series (no need to specify `alphabet_size`)
ts = [1, 2, 1, 2, 1, 1, 1, 2]
compression_complexity(x, EffortToCompress())
```

Multivariate time series (given as [`Dataset`](@ref)s) also work.

```julia
# A multivariate time series from a three-letter alphabet
x1 = rand(0:2, 100)
x2 = rand(0:2, 100)
x3 = rand(0:2, 100)
ts = Dataset(x1, x2, x3)
compression_complexity(ts, EffortToCompress(alphabet_size = 3))
```

Sliding window estimators automatically handle window creation,
and returns a vector of complexity measures computed on each window.

```julia
ts = rand(0:1, 100)
alg = EffortToCompressSlidingWindow(window_size = 10, step = 5)
compression_complexity(x, alg)
```

## Joint 

The joint compression complexity between two time series can also
be computed, but then we must specify the alphabet size to ensure
correctness of the result. Below

```julia
x = rand(0:2, 100)
y = rand(0:2, 100)
alg = EffortToCompress(normalize = true, alphabet_size = 2)
compression_complexity(x, y, alg)
```

This also works for multivariate time series.

```julia
# Variables in the first dataset have an alphabet size of 2
x1, x2 = rand(0:1, 1000), rand(0:1, 1000)

# Variables in the second dataset have an alphabet size of 4
y1, y2, y3 = rand(0:3, 1000), rand(0:3, 1000), rand(0:3, 1000)

d1 = Dataset(x1, x2)
d2 = Dataset(y1, y2, y3)
alg = EffortToCompress(normalize = true)
compression_complexity(d1, d2, alg, 2, 4)
```

See also: [`EffortToCompress`](@ref), [`EffortToCompressSlidingWindow`](@ref).

[^Kathpalia2019]: Kathpalia, A., & Nagaraj, N. (2019). Data-based intervention approach for Complexity-Causality measure. PeerJ Computer Science, 5, e196.
"""
function compression_complexity end
