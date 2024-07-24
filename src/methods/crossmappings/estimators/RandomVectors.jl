using Random

export RandomVectors

"""
    RandomVectors <: CrossmapEstimator
    RandomVectors(definition::CrossmapMeasure; libsizes, replace = false, 
        rng = Random.default_rng())

Cross map *once* over  `N = length(libsizes)` different "point libraries", where 
point indices are selected randomly (not considering time ordering). 

This is method 3 from [Luo2015](@citet). See [`CrossmapEstimator`](@ref) for an in-depth 
explanation of what "library" means in this context.

## Description

The cardinality of the point libraries are given by `libsizes`. One set of 
random point indices is selected per `L ∈ libsizes`, and the `i`-th 
library has cardinality `k = libsizes[i]`. 

Point indices within each library are randomly selected, independently of other libraries.
A user-specified `rng` may be specified for reproducibility. The `replace` argument
controls whether sampling is done with or without replacement. If the time series
you're cross mapping between have length `M`, and `Lᵢ < M` for any `Lᵢ ∈ libsizes`,
then you must set `replace = true`.

## Returns

The return type when used with [`association`](@ref) depends on the type of `libsizes`.
- If `libsizes` is an `Int` (a single library), then a single cross-map estimate is returned.
- If `libsizes` is an `AbstractVector{Int}` (multiple libraries), then a vector of cross-map
    estimates is returned --- one per library.

See also: [`CrossmapEstimator`](@ref).
"""
struct RandomVectors{M <: CrossmapMeasure, I, R} <: CrossmapEstimator{M, I, R}
    definition::M
    libsizes::I
    rng::R
    replace::Bool
    function RandomVectors(definition::M; libsizes::I, replace::Bool = false,
            rng::R = Random.default_rng()) where {M<:CrossmapMeasure, I, R}
        new{M, I, R}(definition, libsizes, rng, replace)
    end
end

function library_indices(est::RandomVectors, i::Int, target, args...)
    N = length(target)
    L = est.libsizes[i]
    return library_indices(measure, est, N, L)
end

function library_indices(est::RandomVectors, N::Int, L::Int)
    return sample(est.rng, 1:N, L; replace = est.replace)
end
