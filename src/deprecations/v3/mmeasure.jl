export m_measure

function m_measure(x, y; kwargs...)
    @warn(
        "Convenience function `m_measure` is deprecated. " *
        "Use `m_measure(MMeasure(; kwargs...), source, target)` instead."
    )
    return association(MMeasure(; kwargs...), x, y)
end

function m_measure(measure::MMeasure, x, y; kwargs...)
    @warn(
        "Convenience function `m_measure` is deprecated. " *
        "Use `association(MMeasure(; kwargs...), source, target) instead."
    )
    return association(MMeasure(; kwargs...), x, y)
end
