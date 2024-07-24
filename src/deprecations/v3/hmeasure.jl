export h_measure

function h_measure(x, y; kwargs...)
    @warn(
        "Convenience function `h_measure` is deprecated. " *
        "Use `h_measure(HMeasure(; kwargs...), source, target)` instead."
    )
    return association(HMeasure(; kwargs...), x, y)
end

function h_measure(measure::HMeasure, x, y; kwargs...)
    @warn(
        "Convenience function `h_measure` is deprecated. " *
        "Use `association(HMeasure(; kwargs...), source, target) instead."
    )
    return association(HMeasure(; kwargs...), x, y)
end
