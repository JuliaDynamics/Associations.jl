export l_measure

function l_measure(x, y; kwargs...)
    @warn(
        "Convenience function `l_measure` is deprecated. " *
        "Use `l_measure(LMeasure(; kwargs...), source, target)` instead."
    )
    return association(LMeasure(; kwargs...), x, y)
end

function l_measure(measure::LMeasure, x, y; kwargs...)
    @warn(
        "Convenience function `l_measure` is deprecated. " *
        "Use `association(LMeasure(; kwargs...), source, target) instead."
    )
    return association(LMeasure(; kwargs...), x, y)
end
