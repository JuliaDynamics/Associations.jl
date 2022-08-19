abstract type CompressionComplexityAlgorithm end

export compression_complexity

"""
    compression_complexity(x, algorithm::EffortToCompress) → N
    compression_complexity(x, algorithm::EffortToCompressSlidingWindow) → Vector{N}

    compression_complexity(x::AbstractVector{Int}, y::AbstractVector{Int}, 
        algorithm::EffortToCompress) → N    
    compression_complexity(x::AbstractVector{Int}, y::AbstractVector{Int}, 
        algorithm::EffortToCompressSlidingWindow) → Vector{N}  

Compute the compression complexity of the pre-binned/pre-symbolized (integer-valued)
univariate (`Vector{Int}`) or multivariate (`Dataset`) time series `x`, using the given 
`algorithm`. 

Or, in the case two symbolized time series `x` and `y` are provided, compute 
the joint compression complexity. For [`EffortToCompress`](@ref), this corresponds to 
the joint ETC described in [^Kathpalia2019].

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

[^Kathpalia2019]: Kathpalia, A., & Nagaraj, N. (2019). Data-based intervention approach for Complexity-Causality measure. PeerJ Computer Science, 5, e196.
"""
function compression_complexity end

include("effort_to_compress/effort_to_compress.jl")


