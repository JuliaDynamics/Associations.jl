using StateSpaceSets: AbstractDataset
using Distances

export DistanceCorrelation
export distance_correlation

"""
    DistanceCorrelation

The distance correlation (Székely et al., 2007)[^Székely2007] measure quantifies
potentially nonlinear associations between pairs of variables.

## Usage

- Use with [`independence`](@ref) to perform a formal hypothesis test for
    pairwise dependence.
- Use with [`distance_correlation`](@ref) to compute the raw distance correlation
    coefficient.
"""
struct DistanceCorrelation <: AssociationMeasure end

"""
    distance_correlation(x, y) → dcor ∈ [0, 1]

Compute the empirical/sample distance correlation (Székely et al., 2007)[^Székely2007],
here called `dcor`, between datasets `x` and `y`.

[^Székely2007]:
    Székely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007). Measuring and testing
    dependence by correlation of distances. The annals of statistics, 35(6), 2769-2794.
"""
function distance_correlation(x::ArrayOrDataset, y::ArrayOrDataset)
    return estimate(DistanceCorrelation(), x, y)
end

# Common interface for higher-level methods.
function estimate(measure::DistanceCorrelation, X, Y)
    # TODO: Future optimization: this could be quicker if we only compute distances once
    # for X and once for Y. Currently, they are computed twice each.
    𝒱ₙ²xy = distance_covariance(X, Y)
    𝒱ₙ²x = distance_covariance(X)
    𝒱ₙ²y = distance_covariance(Y)

    if 𝒱ₙ²x * 𝒱ₙ²y > 0
        return sqrt(𝒱ₙ²xy / sqrt(𝒱ₙ²x * 𝒱ₙ²y))
    else
        return 0.0
    end
end


"""
    distance_covariance(x, y) → dcov::Real

Compute the empirical/sample distance covariance (Székely et al., 2007)[^Székely2007]
between datasets `x` and `y`.

[^Székely2007]:
    Székely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007). Measuring and testing
    dependence by correlation of distances. The annals of statistics, 35(6), 2769-2794.
"""
function distance_covariance(X::ArrayOrDataset, Y::ArrayOrDataset)
    x = Dataset(X)
    y = Dataset(Y)
    Lx, Ly = length(x), length(y)
    Lx == Ly || throw(ArgumentError("Inputs `x` and `y` must have same length"))
    N = length(x)
    # The subscript notation in the paper is a bit messy, but it really just refers
    # to column-wise (āₖs), row-wise (āₗs) and overall (ā) means of a pairwise distance
    # matrix (and analogously for b̄ₖs, b̄ₗs and b̄)
    A = pairwise(Euclidean(), x)
    B = pairwise(Euclidean(), y)
    āₖs = mean(A, dims = 1) # col-wise mean
    āₗs = mean(A, dims = 2) # row-wise mean
    ā = mean(A)
    b̄ₖs = mean(B, dims = 1) # col-wise mean
    b̄ₗs = mean(B, dims = 2) # row-wise mean
    b̄ = mean(B)

    𝒱ₙ² = 0.0
    for l = 1:N
        āₖ = āₖs[l]
        b̄ₖ = b̄ₖs[l]
        for k = 1:N
            Aₖₗ = A[k, l] - āₖ - āₗs[k] + ā
            Bₖₗ = B[k, l] - b̄ₖ - b̄ₗs[k] + b̄
            𝒱ₙ² += Aₖₗ * Bₖₗ
        end
    end
    𝒱ₙ² /= N^2

    return 𝒱ₙ²
end
distance_covariance(x::ArrayOrDataset) = distance_variance(Dataset(x))

"""
    distance_variance(x) → dvar::Real

Compute the empirical/sample distance variance (Székely et al., 2007)[^Székely2007]
for dataset `x`.

[^Székely2007]:
    Székely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007). Measuring and testing
    dependence by correlation of distances. The annals of statistics, 35(6), 2769-2794.
"""
function distance_variance(X::ArrayOrDataset)
    x = Dataset(X)
    N = length(x)
    A = pairwise(Euclidean(), Dataset(x))
    āₖs = mean(A, dims = 1) # col-wise mean
    āₗs = mean(A, dims = 2) # row-wise mean
    ā = mean(A)
    𝒱ₙ² = 0.0
    for l = 1:N
        for k = 1:N
            Aₖₗ = A[k, l] - āₖs[l] - āₗs[k] + ā
            𝒱ₙ² += Aₖₗ^2
        end
    end
    𝒱ₙ² /= N^2

    return 𝒱ₙ²
end
