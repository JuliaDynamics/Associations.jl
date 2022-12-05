"""
    ExpandingSegment <: CrossmapEstimator
    ExpandingSegment(L::Int)

Indicatates that cross mapping is performed on a contiguous time series segment/window,
starting from the first available data point up to the `L`th data point.
"""
struct ExpandingSegment <: CrossmapEstimator
    L::Int
end

ExpandingSegment() = error(crossmap_err())
ExpandingSegment(measure::CrossmapMeasure, x::AbstractVector) =
    ExpandingSegment(max_segmentlength(measure, x))


function predict(measure::CrossmapMeasure, est::ExpandingSegment,
        t::AbstractVector, s::AbstractVector,
        L::Int = max_segmentlength(s, measure))
    # Embed together, so that time indices are aligned properly.
    mixed_embedding = embed_for_crossmap(measure, s, t)
    S̄ = mixed_embedding[1:L, 1:end-1]
    t = mixed_embedding[1:L, end]
    t̂ = predict(measure, S̄[1:L], t[1:L])
    return t, t̂
end
