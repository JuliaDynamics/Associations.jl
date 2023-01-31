export PearsonCorrelation
export pearson_correlation

"""
    PearsonCorrelation

The Pearson correlation of two variables.

This type exists to be used with [`independence`](@ref) testing. If you only need the
correlation coefficient, do [`pearson_correlation`](@ref)`(x, y)`.

## Description

The sample [Pearson correlation coefficient](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient)
for real-valued random variables ``X`` and ``Y`` with associated samples
``\\{x_i\\}_{i=1}^N`` and ``\\{y_i\\}_{i=1}^N`` is defined as

```math
\\rho_{xy} = \\dfrac{\\sum_{i=1}^n (x_i - \\bar{x})(y_i - \\bar{y}) }{\\sqrt{\\sum_{i=1}^N (x_i - \\bar{x})^2}\\sqrt{\\sum_{i=1}^N (y_i - \\bar{y})^2}},
```

where ``\\bar{x}`` and ``\\bar{y}`` are the means of the observations ``x_k`` and ``y_k``,
respectively.
"""
struct PearsonCorrelation end

"""
    pearson_correlation(x::VectorOrDataset, y::VectorOrDataset)

Compute the [`PearsonCorrelation`](@ref) between `x` and `y`, which must each be
1-dimensional.
"""
function pearson_correlation(x::VectorOrDataset, y::VectorOrDataset)
    return estimate(PearsonCorrelation(), x, y)
end

# Common interface for higher-level methods.
function estimate(measure::PearsonCorrelation,
        x::VectorOrDataset{1, T},
        y::VectorOrDataset{1, T}) where T

    Lx, Ly = length(x), length(y)
    Lx == Ly || throw(ArgumentError("Inputs `x` and `y` must have same length"))
    x̄ = extract_mean(x)
    ȳ = extract_mean(y)
    num = 0.0
    den = 0.0
    for (xᵢ, yᵢ) in zip(pt_generator(x), pt_generator(y))
        num += (xᵢ - x̄)*(yᵢ - ȳ)
        den += sqrt((xᵢ - x̄)^2) * sqrt((yᵢ - ȳ)^2)
    end
    ρ = num / den
    return ρ
end

# Silly, but 1-dimensional datasets needs special indexing (because each point is a vector,
# not a value).
pt_generator(x::AbstractDataset{1}) = (x[1] for x in x)
pt_generator(x::AbstractVector) = (x for x in x)

# Ensures type stability.
extract_mean(x::AbstractVector) = mean(x)
extract_mean(x::AbstractDataset{1}) = mean(x)[1]
