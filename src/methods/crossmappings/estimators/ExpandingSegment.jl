export ExpandingSegment

"""
    ExpandingSegment <: CrossmapEstimator
    ExpandingSegment(; libsizes::Int, rng = Random.default_rng())

Indicates that cross mapping is performed on a contiguous time series segment/window,
starting from the first available data point up to the `L`th data point.

If used in an ensemble setting, the estimator is applied to time indices `Lmin:step:Lmax`
of the joint embedding.
"""
struct ExpandingSegment{I, R} <: CrossmapEstimator{I, R}
    libsizes::I
    # For other estimators, `rng` is used for ensemble analyses. For `ExpandingSegment`,
    # an ensemble doesn't make sense, because there is no random sampling involved.
    # However, when the input is uncertain data, the `rng` *does* matter, so we have
    # it here for convenience.
    rng::R

    function ExpandingSegment(; libsizes::I, rng::R = Random.default_rng()) where {I, R}
        new{I, R}(libsizes, rng)
    end
end
