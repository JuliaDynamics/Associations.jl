export pearson_correlation

function pearson_correlation(x, y; kwargs...)
    @warn(
        "Convenience function `pearson_correlation` is deprecated. " *
        "Use `association(PearsonCorrelation(; kwargs...), source, target)` instead."
    )
    association(PearsonCorrelation(; kwargs...), x, y)
end