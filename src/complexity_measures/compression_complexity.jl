abstract type CompressionComplexityAlgorithm end

export compression_complexity

"""
    compression_complexity(x, algorithm::EffortToCompress) → N    
    compression_complexity(x, algorithm::EffortToCompressSlidingWindow) → Vector{N}

Compute the compression complexity of the pre-binned/pre-symbolized (integer-valued)
univariate (`Vector{Int}`) or multivariate (`Dataset`) time series, using the given 
`algorithm`.

See also: [`EffortToCompress`](@ref), [`EffortToCompressSlidingWindow`](@ref).

## Returns

See the documentation for individual algorithms for details on return values.

Some `algorithm`s have a sliding window version that may be optimized for repeated 
application. When using these algorithms, a vector of compression complexity measures - 
one for each window - is returned.

## Examples

### Effort-to-compress

```julia
# A univariate time series (no need to specify `alphabet_size`)
ts = [1, 2, 1, 2, 1, 1, 1, 2]
compression_complexity(x, EffortToCompress())
```

```julia
# A multivariate time series from a three-letter alphabet
ts = rand(0:2, 30)
compression_complexity(x, EffortToCompress(alphabet_size = 3))
```

```julia
ts = rand(0:1, 100)
alg = EffortToCompressSlidingWindow(window_size = 10, step = 5)
compression_complexity(x, alg)
```
"""
function compression_complexity end

include("effort_to_compress/effort_to_compress.jl")


