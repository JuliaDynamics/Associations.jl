using Random

export RandomSegment

"""
    RandomSegment <: CrossmapEstimator
    RandomSegment(definition::CrossmapMeasure; libsizes::Int, rng = Random.default_rng())

Cross map *once* over `N = length(libsizes)` different "point libraries", where 
point indices are selected as time-contiguous segments with random starting points.

This is method 2 from [Luo2015](@cite). See [`CrossmapEstimator`](@ref) for an in-depth 
explanation of what "library" means in this context.

## Description

The cardinality of the point index segments are given by `libsizes`. One segment 
with a randomly selected starting point is picked per `L ∈ libsizes`, and the `i`-th 
point index segment has cardinality `k = libsizes[i]`. 

The starting point for each library is selected independently of other libraries.
A user-specified `rng` may be specified for reproducibility. If the time series
you're cross mapping between have length `M`, and `Lᵢ < M` for any `Lᵢ ∈ libsizes`,
then an error will be thrown.

A user-specified `rng` may be specified for reproducibility.

## Returns

The return type when used with [`association`](@ref) depends on the type of `libsizes`.
- If `libsizes` is an `Int` (a single library), then a single cross-map estimate is returned.
- If `libsizes` is an `AbstractVector{Int}` (multiple libraries), then a vector of cross-map
    estimates is returned --- one per library.

See also: [`CrossmapEstimator`](@ref).
"""
struct RandomSegment{M <: CrossmapMeasure, I, R} <: CrossmapEstimator{M, I, R}
    definition::M
    libsizes::I
    rng::R
    function RandomSegment(definition::M; libsizes::I, rng::R = Random.default_rng()) where {M <: CrossmapMeasure, I, R}
        new{M, I, R}(definition, libsizes, rng)
    end
end

function library_indices(est::RandomSegment, i::Int, target, args...)
    definition = est.definition
    N = length(target)
    L = est.libsizes[i]
    Lmax = max_segmentlength(definition, target)
    L <= Lmax ||
        throw(ArgumentError("L ($L) > Lmax ($Lmax). Use a smaller segment length (some points are lost when embedding)."))
    return library_indices(definition, est, N, L)
end

function library_indices(est::RandomSegment, N::Int, L::Int)
    startidx = sample(est.rng, 1:(N - L)) # random segment starting point
    return startidx:startidx+L-1
end