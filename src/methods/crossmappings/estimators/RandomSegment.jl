using Random

export RandomSegment

"""
    RandomSegment <: CrossmapEstimator
    RandomSegment(; libsizes::Int, rng = Random.default_rng())

Indicatates that cross mapping is performed on contiguous time series
segments/windows of length `L` with a randomly selected starting point.

This is method 2 from [Luo2015](@cite).
"""
struct RandomSegment{I, R} <: CrossmapEstimator{I, R}
    libsizes::I
    rng::R
    function RandomSegment(; libsizes::I, rng::R = Random.default_rng()) where {I, R}
        new{I, R}(libsizes, rng)
    end
end
