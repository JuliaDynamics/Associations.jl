export distance_correlation

function distance_correlation(x, y; kwargs...)
    @warn(
        "Convenience function `distance_correlation` is deprecated. " *
        "Use `association(DistanceCorrelation(), x, y)` instead."
    )
    association(DistanceCorrelation(; kwargs...), x, y)
end

function distance_correlation(x, y, z; kwargs...)
    @warn(
        "Convenience function `distance_correlation` is deprecated. " *
        "Use `association(DistanceCorrelation(), x, y, z)` instead."
    )
    association(DistanceCorrelation(; kwargs...), x, y, z)
end
