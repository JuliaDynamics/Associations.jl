using Random

export RandomSegment

"""
    RandomSegment <: CrossmapEstimator
    RandomSegment(definition::CrossmapMeasure; 
        libsizes::Int, rng = Random.default_rng())

Indicatates that cross mapping is performed on contiguous time series
segments/windows of length `L` with a randomly selected starting point.

This is method 2 from [Luo2015](@cite).
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