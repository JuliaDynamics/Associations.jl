export PearsonCorrelation

"""
    PearsonCorrelation

The Pearson correlation of two variables.

## Usage

- Use with [`association`](@ref) to compute the raw Pearson correlation coefficient.
- Use with [`independence`](@ref) to perform a formal hypothesis test for pairwise dependence
    using the Pearson correlation coefficient.

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
struct PearsonCorrelation <: CorrelationMeasure end

# Common interface for higher-level methods.
function association(measure::PearsonCorrelation,
        x::VectorOrStateSpaceSet{1, T},
        y::VectorOrStateSpaceSet{1, T}) where T
    Lx, Ly = length(x), length(y)
    Lx == Ly || throw(ArgumentError("Inputs `x` and `y` must have same length"))
    x̄ = extract_mean(x)
    ȳ = extract_mean(y)
    num = 0.0
    for (xᵢ, yᵢ) in zip(pt_generator(x), pt_generator(y))
        num += (xᵢ - x̄)*(yᵢ - ȳ)
    end
    den = sqrt(
        sum((xᵢ - x̄)^2 for xᵢ in pt_generator(x)) *
        sum((yᵢ - ȳ)^2 for yᵢ in pt_generator(y))
    )
    ρ = num / den
    return ρ
end

function association(measure::PearsonCorrelation, est::Nothing, x, y)
    return association(measure, x, y)
end

# Silly, but 1-dimensional StateSpaceSets needs special indexing (because each point is a vector,
# not a value).
pt_generator(x::AbstractStateSpaceSet{1}) = (x[1] for x in x)
pt_generator(x::AbstractVector) = (x for x in x)

# Ensures type stability.
extract_mean(x::AbstractVector) = mean(x)
extract_mean(x::AbstractStateSpaceSet{1}) = mean(x)[1]
