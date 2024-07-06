using Neighborhood: searchstructure, bulkisearch
using Neighborhood: NeighborNumber, Theiler
using StateSpaceSets: AbstractStateSpaceSet
using Distances: SqEuclidean, Euclidean
using Distances: pairwise, evaluate

export m_measure
export MMeasure

"""
    MMeasure <: AssociationMeasure
    MMeasure(; K::Int = 2, dx = 2, dy = 2, τx = - 1, τy = -1, w = 0)

The `MMeasure` [Andrzejak2003](@cite) is a pairwise association
measure. It quantifies the probability with which close state of a target
timeseries/embedding are mapped to close states of a source timeseries/embedding.

Note that `τx` and `τy` are negative by convention. See docstring for [`SMeasure`](@ref)
for an explanation.

## Usage

- Use with [`association`](@ref)/[`m_measure`](@ref) to compute the raw m-measure statistic.
- Use with [`independence`](@ref) to perform a formal hypothesis test for directional dependence.

## Description

The `MMeasure` is based on [`SMeasure`](@ref) and [`HMeasure`](@ref). It is given by

```math
M^{(k)}(x|y) = \\dfrac{1}{N} \\sum_{i=1}^{N}
\\log \\left( \\dfrac{R_i(x) - R_i^{(k)}(x|y)}{R_i(x) - R_i^k(x)} \\right),
```

where ``R_i(x)`` is computed as for [`HMeasure`](@ref), while ``R_i^k(x)`` and
``R_i^{(k)}(x|y)`` is computed as for [`SMeasure`](@ref).
Parameters also have the same meaning as for [`SMeasure`](@ref)/[`HMeasure`](@ref).
"""
Base.@kwdef struct MMeasure{M, TM} <: AssociationMeasure
    K::Int = 2
    metric::M = SqEuclidean()
    tree_metric::TM = Euclidean()
    τx::Int = -1
    τy::Int = -1
    dx::Int = 2
    dy::Int = 2
    w::Int = 0
end

"""
    m_measure(measure::MMeasure, x::VectorOrStateSpaceSet, y::VectorOrStateSpaceSet)

Compute the [`MMeasure`](@ref) from source `x` to target `y`.
"""
function m_measure(measure::MMeasure, x::VectorOrStateSpaceSet, y::VectorOrStateSpaceSet)
    return association(measure, x, y)
end

# Internal method for use with `independence`
function association(measure::MMeasure, x::AbstractStateSpaceSet, y::AbstractStateSpaceSet)
    (; K, metric, tree_metric, τx, τy, dx, dy, w) = measure

    # Match length of StateSpaceSets by excluding end points.
    lx = length(x); ly = length(y)
    lx > ly ? X = x[1:ly, :] : X = x
    ly > lx ? Y = y[1:lx, :] : Y = y
    N = length(X)

    T = eltype(1.0)
     # Pre-allocate vectors to hold indices and distances during loops
    dists_x = zeros(T, K)
    dists_x_cond_y = zeros(T, K)

    # Rᵢx := mean squared distance to all other points
    # Rᵢᵏx := Mean squared distances in X, and
    # Rᵢᵏxy := mean squared distances in X conditioned on Y
    Rᵢx = zeros(T, N)
    Rᵢᵏx = zeros(T, N)
    Rᵢᵏxy = zeros(T, N)

    # Search for the K nearest neighbors of each points in both X and Y
    treeX = searchstructure(KDTree, X, tree_metric)
    treeY = searchstructure(KDTree, Y, tree_metric)
    dx = pairwise(metric, X)
    neighborhoodtype, theiler = NeighborNumber(K), Theiler(w)
    idxs_X = bulkisearch(treeX, X, neighborhoodtype, theiler)
    idxs_Y = bulkisearch(treeY, Y, neighborhoodtype, theiler)

    for n in 1:N
        pxₙ = X[n]

        for j = 1:K
            rₙⱼ = idxs_X[n][j] # nearest neighbor indices in X
            sₙⱼ = idxs_Y[n][j] # nearest neighbor indices in Y
            dists_x[j] = evaluate(metric, pxₙ, X[rₙⱼ])
            dists_x_cond_y[j] = evaluate(metric, pxₙ, X[sₙⱼ])
        end
        Rᵢx[n] = sum(dx[:, n]) / N
        Rᵢᵏx[n] = sum(dists_x) / K
        Rᵢᵏxy[n] = sum(dists_x_cond_y) / K
    end

    return sum((Rᵢx .- Rᵢᵏxy) ./ (Rᵢx .- Rᵢᵏx)) / N
end
