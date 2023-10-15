# Full names
@deprecate PMI PartialMutualInformation

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
    msg = "crossmap(measure::CrossmapMeasure, est::CrossmapEstimator, args...)" * 
        " is no longer possible. The signature is now " * 
        "crossmap(est::CrossmapEstimator{<:CrossmapMeasure}, args...), so you must give the " * 
        "measure definition as the first argument to the estimator." 
    throw(ArgumentError(msg))
end

function crossmap(x::AbstractVector, y::AbstractVector, args...)
    msg = "crossmap(x::AbstractVector, y::AbstractVector, args...)" * 
    " is no longer possible. The signature is now " * 
    "crossmap(est::CrossmapEstimator{<:CrossmapMeasure}, args...). " *
    "For example, you can do `crossmap(RandomSegment(CCM(); libsizes = 10:10:50), x, y)`" 
    throw(ArgumentError(msg))
end

function pai(x::AbstractVector, y::AbstractVector, args...)
    msg = "pai(x::AbstractVector, y::AbstractVector, args...)" * 
    " is no longer possible. The signature is now " * 
    "pai(est::CrossmapEstimator{<:CrossmapMeasure}, args...). " *
    "For example, you can do `pai(RandomSegment(CCM(); libsizes = 10:10:50), x, y)`" 
    throw(ArgumentError(msg))
end


    # function jdd(source, target; kw...)
    #     if !isempty(kw)
    #         @warn(
    #             "Providing keywords to `jdd` is deprecated. " *
    #             "Use `jdd(JointDistanceDistribution(; kwargs...), source, target) instead of " *
    #             "`jdd(source, target; kwargs...)`"
    #         )
    #     end
    #     estimate(JointDistanceDistribution(; kw...), source, target)
    # end