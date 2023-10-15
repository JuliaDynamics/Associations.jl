export ExpandingSegment

"""
    ExpandingSegment <: CrossmapEstimator
    ExpandingSegment(definition::CrossmapMeasure; 
        libsizes::Int, rng = Random.default_rng())

Indicatates that cross mapping is performed on a contiguous time series segment/window,
starting from the first available data point up to the `L`th data point.

If used in an ensemble setting, the estimator is applied to time indices `Lmin:step:Lmax`
of the joint embedding.
"""
struct ExpandingSegment{M <: CrossmapMeasure, I, R} <: CrossmapEstimator{M, I, R}
    definition::M
    libsizes::I
    # For other estimators, `rng` is used for ensemble analyses. For `ExpandingSegment`,
    # an ensemble doesn't make sense, because there is no random sampling involved.
    # However, when the input is uncertain data, the `rng` *does* matter, so we have
    # it here for convenience.
    rng::R

    function ExpandingSegment(definition::M; libsizes::I, rng::R = Random.default_rng()) where {M <: CrossmapMeasure, I, R}
        new{M, I, R}(definition, libsizes, rng)
    end
end

function library_indices(est::ExpandingSegment, i::Int, target, args...)
    definition = est.definition
    Lmax = max_segmentlength(est.definition, target)
    L = est.libsizes[i]
    L <= Lmax || throw(ArgumentError("L ($L) > Lmax ($Lmax). Use a smaller segment length (some points are lost when embedding)."))
    return library_indices(est, length(target), L)
end
function library_indices(est::ExpandingSegment, N::Int, L::Int)
    return 1:L
end
