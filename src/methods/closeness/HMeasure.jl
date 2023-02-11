using NearestNeighbors, Distances, Neighborhood, DelayEmbeddings
using Distances: SqEuclidean, Euclidean
export h_measure
export HMeasure

"""
    HMeasure <: AssociationMeasure
    HMeasure(; K::Int = 2, dx::Int = 2, my::Int = 2, τx::Int = 1, τy::Int = 1, w::Int = 0)

The `HMeasure` (Grassberger et al., 1999)[^Grassberger1999] is a pairwise association
measure. It quantifies the probability with which close state of a target
timeseries/embedding are mapped to close states of a source timeseries/embedding.

## Usage

- Use with [`independence`](@ref) to perform a formal hypothesis test for directional dependence.
- Use with [`h_measure`](@ref) to compute the raw h-measure statistic.

## Description

The `HMeasure` (Grassberger et al., 1999)[^Grassberger1999] is similar to the
[`SMeasure`](@ref), but the numerator of the formula is replaced by ``R_i(x)``, the mean
squared Euclidean distance to *all other points*, and there is a ``\\log``-term inside
the sum:

```math
H^{(k)}(x|y) = \\dfrac{1}{N} \\sum_{i=1}^{N}
\\log \\left( \\dfrac{R_i(x)}{R_i^{(k)}(x|y)} \\right).
```

Parameters are the same and ``R_i^{(k)}(x|y)`` is computed as for [`SMeasure`](@ref).

[^Grassberger1999]:
    Arnhold, J., Grassberger, P., Lehnertz, K., & Elger, C. E. (1999). A robust method for detecting
    interdependences: application to intracranially recorded EEG. Physica D:
    Nonlinear Phenomena, 134(4), 419-430.
"""
Base.@kwdef struct HMeasure{M, TM} <: AssociationMeasure
    K::Int = 2
    metric::M = SqEuclidean()
    tree_metric::TM = Euclidean()
    τx::Int = 1
    τy::Int = 1
    dx::Int = 2
    dy::Int = 2
    w::Int = 0
end

"""
    h_measure(measure::HMeasure, x::VectorOrDataset, y::VectorOrDataset)

Compute the [`HMeasure`](@ref) from source `x` to target `y`.
"""
function h_measure(measure::HMeasure, x::VectorOrDataset, y::VectorOrDataset)
    return estimate(measure, x, y)
end

# Internal method for use with `independence`
function estimate(measure::HMeasure, x::AbstractDataset, y::AbstractDataset)
    (; K, metric, tree_metric, τx, τy, dx, dy, w) = measure

    # Match length of datasets by excluding end points.
    lx = length(x); ly = length(y)
    lx > ly ? X = x[1:ly, :] : X = x
    ly > lx ? Y = y[1:lx, :] : Y = y
    N = length(X)

    T = eltype(1.0)
     # Pre-allocate vectors to hold indices and distances during loops
    dists_x_cond_y = zeros(T, K)

    # Rᵢx := Mean squared distance to all other points in X, and
    # Rᵢᵏxy := mean squared distances in X conditioned on Y
    Rᵢx = zeros(T, N)
    Rᵢᵏxy = zeros(T, N)

    # Search for the K nearest neighbors of each points in both X and Y
    dx = pairwise(metric, X)
    treeY = searchstructure(KDTree, Y, tree_metric)
    neighborhoodtype, theiler = NeighborNumber(K), Theiler(w)
    idxs_Y = bulkisearch(treeY, Y, neighborhoodtype, theiler)

    for n in 1:N
        pxₙ = X[n]
        for j = 1:K
            sₙⱼ = idxs_Y[n][j] # nearest neighbor indices in Y
            dists_x_cond_y[j] = evaluate(metric, pxₙ, X[sₙⱼ])
        end
        Rᵢᵏxy[n] = sum(dists_x_cond_y) / K
        Rᵢx[n] = sum(dx[:, n]) / N
    end

    return sum(log.(Rᵢx ./ Rᵢᵏxy)) / N
end
