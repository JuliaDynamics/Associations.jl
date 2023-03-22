using RecurrenceAnalysis: RecurrenceMatrix
export RMCD
export rmcd

"""
    RMCD <: AssociationMeasure
    RMCD(; r, metric = Euclidean(), base = 2)

The recurrence measure of conditional dependence, or RMCD (Ramos et al., 2017)[^Ramos2017],
is a recurrence-based measure that mimics the conditional mutual
information, but uses recurrence probabilities.

`r` is a mandatory keyword which specifies the recurrence threshold when constructing
recurrence matrices. It can be instance of
any subtype of `AbstractRecurrenceType` from
[RecurrenceAnalysis.jl](https://juliadynamics.github.io/RecurrenceAnalysis.jl/stable/).
To use any `r` that is not a real number, you have to do `using RecurrenceAnalysis` first.
The `metric` is any valid metric
from [Distances.jl](https://github.com/JuliaStats/Distances.jl).

Both the pairwise and conditional RMCD is non-negative, but due to round-off error,
negative values may occur. If that happens, an RMCD value of `0.0` is returned.

## Usage

- Use with [`independence`](@ref) to perform a formal hypothesis test for pairwise
    or conditional association.
- Use with [`rmcd`](@ref) to compute the raw RMCD for pairwise or conditional association.

## Description

The RMCD measure is defined by

```math
I_{RMCD}(X; Y | Z) = \\dfrac{1}{N}
\\sum_{i} \\left[
\\dfrac{1}{N} \\sum_{j} R_{ij}^{X, Y, Z}
\\log \\left(
    \\dfrac{\\sum_{j} R_{ij}^{X, Y, Z} \\sum_{j} R_{ij}^{Z} }{\\sum_{j} \\sum_{j} R_{ij}^{X, Z} \\sum_{j} \\sum_{j} R_{ij}^{Y, Z}}
    \\right)
\\right],
```

where  `base` controls the base of the logarithm.
``I_{RMCD}(X; Y | Z)`` is zero when ``Z = X``, ``Z = Y`` or
when ``X``, ``Y`` and ``Z`` are mutually independent.

Our implementation allows dropping the third/last argument, in which
case the following mutual information-like quantitity is computed (not
discussed in Ramos et al., 2017).

```math

I_{RMCD}(X; Y) = \\dfrac{1}{N}
\\sum_{i} \\left[
\\dfrac{1}{N} \\sum_{j} R_{ij}^{X, Y}
\\log \\left(
    \\dfrac{\\sum_{j} R_{ij}^{X}  R_{ij}^{Y} }{\\sum_{j} R_{ij}^{X, Y}}
    \\right)
\\right]
```

[^Ramos2017]:
    Ramos, A. M., Builes-Jaramillo, A., Poveda, G., Goswami, B., Macau, E. E.,
    Kurths, J., & Marwan, N. (2017). Recurrence measure of conditional dependence and
    applications. Physical Review E, 95(5), 052206.
"""
Base.@kwdef struct RMCD{R, M, B} <: AssociationMeasure
    r::R
    metric::M = Euclidean()
    base::B = 2
end

"""
    rmcd(measure::RMCD, x, y)
    rmcd(measure::RMCD, x, y, [z, ...])

Estimate the recurrence-based `measure` of dependence between
`x` and `y`, conditional on `z` if given.

Parameters for recurrence matrix estimation are given as a [`RMCD`](@ref) instance.
Inputs `x`, `y`, `z` can be either univariate timeseries or multivariate
[`StateSpaceSet`](@ref)s.
"""
rmcd(measure::RMCD, args...) = estimate(measure, args...)

# For compatibility with independence testing framework.
function estimate(measure::RMCD, est::Nothing, x, y, z)
    return estimate(measure, x, y, z)
end
function estimate(measure::RMCD, est::Nothing, x, y)
    return estimate(measure, x, y)
end

function estimate(measure::RMCD, x, y, z)
    (; r, metric, base) = measure
    @assert length(x) == length(y) == length(z)
    N = length(x)
    rX = RecurrenceMatrix(x, r; metric)
    rY = RecurrenceMatrix(y, r; metric)
    rZ = RecurrenceMatrix(z, r; metric)
    rXZ = rX .* rZ
    rYZ = rY .* rZ
    rXYZ = rX .* rY .* rZ

    Irmcd = 0.0
    # The recurrence matrices are not symmetric for all recurrence types, so we stick
    # to the original indexing in the paper. This could have been optimized by flipping
    # dimensions if symmetric matrices were guaranteed
    for i = 1:N
        joint = @views sum(rXYZ[i, :])
        marginals = @views (joint * sum(rZ[i, :])) / (sum(rXZ[i, :]) * sum(rYZ[i, :]))
        if marginals != 0

        for j = 1:N
                Irmcd += rXYZ[i, j] * log(base, marginals)
            end
        end
        Irmcd /= N
    end
    Irmcd /= N
    # Handle round-off error
    if Irmcd < 0.0
        return 0.0
    else
        return Irmcd
    end
end

# Similar, but analogous to mutual information
function estimate(measure::RMCD, x, y)
    (; r, metric, base) = measure
    @assert length(x) == length(y)
    N = length(x)
    rX = RecurrenceMatrix(x, r; metric)
    rY = RecurrenceMatrix(y, r; metric)
    rXY = rX .* rY

    Irmd = 0.0
    # The recurrence matrices are not symmetric for all recurrence types, so we stick
    # to the original indexing in the paper. This could have been optimized by flipping
    # dimensions if symmetric matrices were guaranteed
    for i = 1:N
        joint = @views sum(rXY[i, :])
        marginals = @views (sum(rX[i, :]) * sum(rY[i, :])) / joint
        for j = 1:N
            if marginals != 0
                Irmd += rXY[i, j] * log(base, marginals)
            end
        end
        Irmd /= N
    end
    Irmd /= N

    # Handle round-off error
    if Irmd < 0.0
        return 0.0
    else
        return Irmd
    end
end
