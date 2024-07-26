export partial_correlation

function partial_correlation(x, y, z; kwargs...)
    @warn(
        "Convenience function `partial_correlation` is deprecated. " *
        "Use `association(PartialCorrelation(), x, y, z)` instead."
    )
    association(PartialCorrelation(; kwargs...), x, y, z)
end