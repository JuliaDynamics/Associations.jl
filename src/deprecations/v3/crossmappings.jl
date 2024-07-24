export pai
export crossmap


# These are not actual deprecations, but breaking changes. Better to error with explicit
# change message.
function RandomVectors(; kwargs...)
    msg = "RandomVectors now takes a `CrossmapMeasure` as the first argument. " *
        "Do `RandomVectors(CCM(); kwargs...)` instead of `RandomVectors(; kwargs...)."
    throw(ArgumentError(msg))
end

function RandomSegment(; kwargs...)
    msg = "RandomSegment now takes a `CrossmapMeasure` as the first argument. " *
        "Do `RandomSegment(CCM(); kwargs...)` instead of `RandomSegment(; kwargs...)."
    throw(ArgumentError(msg))
end

function ExpandingSegment(; kwargs...)
    msg = "ExpandingSegment now takes a `CrossmapMeasure` as the first argument. " *
        "Do `ExpandingSegment(CCM(); kwargs...)` instead of `ExpandingSegment(; kwargs...)."
    throw(ArgumentError(msg))
end

function crossmap(measure::CCMLike, est::CrossmapEstimator, args...)
    msg = "crossmap(measure::CrossmapMeasure, est::CrossmapEstimator, args...) is deprecated. " * 
        "Use `association(est::CrossmapEstiamator, x, y)` instead.`" 
    throw(ArgumentError(msg))
end

function crossmap(x::AbstractVector, y::AbstractVector, args...)
    msg = "crossmap(x::AbstractVector, y::AbstractVector, args...) is deprecated" * 
    "Use `association(RandomSegment(CCM(); libsizes = 10:10:50), x, y)` instead.`" 
    throw(ArgumentError(msg))
end

function pai(x::AbstractVector, y::AbstractVector, args...)
    msg = "pai(x::AbstractVector, y::AbstractVector, args...) is deprecated. " * 
    "Use `association(RandomSegment(CCM(); libsizes = 10:5:30), x, y)` instead.`" 
    throw(ArgumentError(msg))
end
