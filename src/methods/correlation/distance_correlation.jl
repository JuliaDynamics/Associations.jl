using StateSpaceSets: AbstractStateSpaceSet
using Distances
using LinearAlgebra

export DistanceCorrelation
export distance_correlation

"""
    DistanceCorrelation

The distance correlation (Székely et al., 2007)[^Székely2007] measure quantifies
potentially nonlinear associations between pairs of variables. If applied to
three variables, the partial distance correlation (Székely and Rizzo, 2014)[^Székely2014]
is computed.

## Usage

- Use with [`independence`](@ref) to perform a formal hypothesis test for
    pairwise dependence.
- Use with [`distance_correlation`](@ref) to compute the raw distance correlation
    coefficient.

!!! warn
    A partial distance correlation `I = distance_correlation(X, Y, Z)` doesn't
    always guarantee conditional independence `X ⫫ Y | Z`. See Székely and Rizzo (2014)
    for in-depth discussion.

[^Székely2007]:
    Székely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007). Measuring and testing
    dependence by correlation of distances. The annals of statistics, 35(6), 2769-2794.
[^Székely2014]:
    Székely, G. J., & Rizzo, M. L. (2014). Partial distance correlation with methods for
    dissimilarities.
"""
struct DistanceCorrelation <: AssociationMeasure end

"""
    distance_correlation(x, y) → dcor ∈ [0, 1]
    distance_correlation(x, y, z) → pdcor ∈ [0, 1]

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

function ucenter(a) # input `a` is a symmetric distance matrix
    N = size(a, 1)
    Ã  = zeros(size(a))
    f = 1 / (N - 2)
    for k = 1:N
        for l = 1:N
            if k != l
                Ã[k, l] = a[k, l] -
                    f*sum(a[:, l]) -
                    f*sum(a[k, :]) +
                    1 / ((N - 1) * (N - 2)) * sum(a)
            end
        end
    end
    return Ã
end

function inner_prod(a, b)
    size(a) == size(b) || throw(ArgumentError("Matrices must have same size."))
    N = size(a, 1)
    ab = 0.0
    for j = 1:N
        for i = 1:N
            if j != i
                ab += a[i, j] * b[i, j]
            end
        end
    end
    return 1 / (N * (N - 3)) * ab
end

function distance_covariance(X::ArrayOrStateSpaceSet, Y::ArrayOrStateSpaceSet, Z::ArrayOrStateSpaceSet)
    Lx, Ly, Lz = length(X), length(Y), length(Z)
    Lx == Ly == Lz || throw(ArgumentError("Input X, Y and Z must have same lengths."))
    N = Lx
    N >= 4 || throw(ArgumentError("Partial distance covariance is defined for 4 or more points. Got $N"))

    Xds = pairwise(Euclidean(), StateSpaceSet(X))
    Yds = pairwise(Euclidean(), StateSpaceSet(Y))
    Zds = pairwise(Euclidean(), StateSpaceSet(Z))
    Ã = ucenter(Xds)
    B̃ = ucenter(Yds)
    C̃ = ucenter(Zds)

    C̃dotC̃ = inner_prod(C̃, C̃)
    if C̃dotC̃ == 0
        PzX = Ã
        PzY = B̃
    else
        # Orthogonal projection of Ã(x) onto C̃(z)^{⊥}
        PzX = Ã - inner_prod(Ã, C̃) / (C̃dotC̃)*C̃
        # Orthogonal projection of B̃(x) onto C̃(z)^{⊥}
        PzY = B̃ - inner_prod(B̃, C̃) / (C̃dotC̃)*C̃
    end

    return inner_prod(PzX, PzY)
end


# Common interface for higher-level methods.
function estimate(measure::DistanceCorrelation, X, Y, Z)
    Lx, Ly, Lz = length(X), length(Y), length(Z)
    Lx == Ly == Lz || throw(ArgumentError("Input X, Y and Z must have same lengths."))
    N = Lx
    N >= 4 || throw(ArgumentError("Partial distance correlation is defined for 4 or more points. Got $N"))

    Xds = pairwise(Euclidean(), StateSpaceSet(X))
    Yds = pairwise(Euclidean(), StateSpaceSet(Y))
    Zds = pairwise(Euclidean(), StateSpaceSet(Z))
    Ã = ucenter(Xds)
    B̃ = ucenter(Yds)
    C̃ = ucenter(Zds)

    C̃dotC̃ = inner_prod(C̃, C̃)
    if C̃dotC̃ ≈ 0
        PzX = Ã
        PzY = B̃
    else
        # Orthogonal projection of Ã(x) onto C̃(z)^{⊥}
        PzX = Ã - inner_prod(Ã, C̃) / (C̃dotC̃)*C̃
        # Orthogonal projection of B̃(x) onto C̃(z)^{⊥}
        PzY = B̃ - inner_prod(B̃, C̃) / (C̃dotC̃)*C̃
    end

    ip = inner_prod(PzX, PzY)
    norm_XontoZ = inner_prod(PzX, PzX)^(0.5)
    norm_YontoZ = inner_prod(PzY, PzY)^(0.5)
    if norm_XontoZ * norm_YontoZ != 0.0
        return ip / (norm_XontoZ * norm_YontoZ)
    else
        return 0.0
    end
end
