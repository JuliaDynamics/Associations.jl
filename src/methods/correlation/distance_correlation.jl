using StateSpaceSets: AbstractStateSpaceSet
using Distances
using LinearAlgebra

export DistanceCorrelation
export distance_correlation

"""
    DistanceCorrelation

The distance correlation (Székely et al., 2007)[^Székely2007] measure quantifies
potentially nonlinear associations between pairs of variables. If applied to
three variables, the partial distance correlation is computed.

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
here called `dcor`, between StateSpaceSets `x` and `y`.

[^Székely2007]:
    Székely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007). Measuring and testing
    dependence by correlation of distances. The annals of statistics, 35(6), 2769-2794.
"""
function distance_correlation(x::ArrayOrStateSpaceSet, y::ArrayOrStateSpaceSet)
    return estimate(DistanceCorrelation(), x, y)
end

function distance_correlation(x::ArrayOrStateSpaceSet, y::ArrayOrStateSpaceSet,
        z::ArrayOrStateSpaceSet)
    return estimate(DistanceCorrelation(), x, y, z)
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
between StateSpaceSets `x` and `y`.

[^Székely2007]:
    Székely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007). Measuring and testing
    dependence by correlation of distances. The annals of statistics, 35(6), 2769-2794.
"""
function distance_covariance(X::ArrayOrStateSpaceSet, Y::ArrayOrStateSpaceSet)
    x = StateSpaceSet(X)
    y = StateSpaceSet(Y)
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
distance_covariance(x::ArrayOrStateSpaceSet) = distance_variance(StateSpaceSet(x))

"""
    distance_variance(x) → dvar::Real

Compute the empirical/sample distance variance (Székely et al., 2007)[^Székely2007]
for StateSpaceSet `x`.

[^Székely2007]:
    Székely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007). Measuring and testing
    dependence by correlation of distances. The annals of statistics, 35(6), 2769-2794.
"""
function distance_variance(X::ArrayOrStateSpaceSet)
    x = StateSpaceSet(X)
    N = length(x)
    A = pairwise(Euclidean(), StateSpaceSet(x))
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

export ucenter

function ucenter(x::ArrayOrStateSpaceSet) # operates on points
    length(x) >= 4 || throw(ArgumentError("Partial distance correlation is defined for 4 or more points. Got $(length(x))"))
    ds = pairwise(Euclidean(), StateSpaceSet(x))
    return ucenter_distancematrix(ds)
end

function ucenter_distancematrix(ds)
    N = size(ds, 1)
    āₖs = mean(ds, dims = 1) # col-wise mean
    āₗs = mean(ds, dims = 2) # row-wise mean
    Aₖₗ = zeros(size(ds))
    ā = mean(ds)
    for l = 1:N
        for k = 1:N
            Aₖₗ[k, l] = ds[k, l] - āₖs[l] - āₗs[k] + ā
        end
    end
    return Aₖₗ
end


function compute_prod(X̃, Ỹ)
    size(X̃) == size(Ỹ) || throw(ArgumentError("Matrices must have same size."))
    N = size(X̃, 1)
    X̃Ỹ = 0.0
    for j = 1:N
        for i = 1:N
            if j != i
                X̃Ỹ += X̃[i, j] * Ỹ[i, j]
            end
        end
    end
    return X̃Ỹ * 1/(N*(N-3))
end

# Common interface for higher-level methods.
function estimate(measure::DistanceCorrelation, X, Y, Z)
    Lx, Ly, Lz = length(X), length(Y), length(Z)
    Lx == Ly == Lz || throw(ArgumentError("Input X, Y and Z must have same lengths."))
    N = Lx

    Ã = ucenter(X)
    B̃ = ucenter(Y)
    C̃ = ucenter(Z)
    ÃB̃ = compute_prod(Ã, B̃)
    ÃC̃ = compute_prod(Ã, C̃)
    B̃C̃ = compute_prod(B̃, C̃)
    Ãsq = sqrt(Ã ⋅ Ã)
    B̃sq = sqrt(B̃ ⋅ B̃)
    C̃sq = sqrt(C̃ ⋅ C̃)

    Rxy = ÃB̃ / (Ãsq ⋅ B̃sq)
    Rxz = ÃC̃ / (Ãsq ⋅ C̃sq)
    Ryz = B̃C̃ / (B̃sq ⋅ C̃sq)

    if Rxy ^ 2 != 1.0 && Ryz ^ 2 != 1.0
        return (Rxy - Rxz * Rxy) / (sqrt(1 - Rxz^2) * sqrt(1 - Ryz^2))
    else
        PzX = Ã - (Ã ⋅ C̃) / (C̃ ⋅ C̃)*C̃
        PzY = B̃ - (B̃ ⋅ C̃) / (C̃ ⋅ C̃)*C̃
        # Assuming notation in paper means determinant.
        PzXsq = det(PzX)
        PzYsq = det(PzY)
        if PzXsq * PzYsq != 0.0
            return (PzX ⋅ PzY) / (PzXsq * PzYsq)
        else
            return 0.0
        end
    end
end
