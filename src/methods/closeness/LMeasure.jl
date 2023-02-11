using NearestNeighbors, Distances, Neighborhood, StateSpaceSets
using Distances: SqEuclidean, Euclidean

export LMeasure
export l_measure

"""
    LMeasure <: AssociationMeasure
    LMeasure(; K::Int = 2, dx::Int = 2, my::Int = 2, τx::Int = 1, τy::Int = 1)

The `LMeasure` (Chicharro & Andrzejak, 2009)[^^Chicharro20093] is a pairwise association
measure. It quantifies the probability with which close state of a target
timeseries/embedding are mapped to close states of a source timeseries/embedding.

## Usage

- Use with [`independence`](@ref) to perform a formal hypothesis test for directional dependence.
- Use with [`l_measure`](@ref) to compute the raw l-measure statistic.

## Description

`LMeasure` is similar to [`MMeasure`](@ref), but uses distance ranks instead of the raw
distances.

Let ``\\bf{x_i}`` be an embedding vector, and let ``g_{i,j}`` denote the rank that
the distance between ``\\bf{x_i}`` and some other vector ``\\bf{x_j}`` in a sorted
ascending list of distances between ``\\bf{x_i}`` and ``\\bf{x_{i \\neq j}}``
In other words, ``g_{i,j}`` this is just the ``N-1`` nearest neighbor distances sorted )

`LMeasure` is then defined as

```math
L^{(k)}(x|y) = \\dfrac{1}{N} \\sum_{i=1}^{N}
\\log \\left( \\dfrac{G_i(x) - G_i^{(k)}(x|y)}{G_i(x) - G_i^k(x)} \\right),
```

where ``G_i(x) = \\frac{N}{2}`` and ``G_i^K(x) = \\frac{k+1}{2}`` are the mean
and minimal rank, respectively.

The ``y``-conditioned mean rank is defined as

```math
G_i^{(k)}(x|y) = \\dfrac{1}{K}\\sum_{j=1}^{K} g_{i,w_{i, j}},
```

where ``w_{i,j}`` is the index of the ``j``-th nearest neighbor of ``\\bf{y_i}``.

[^Chicharro2009]:
    Chicharro, D., & Andrzejak, R. G. (2009). Reliable detection of directional couplings
    using rank statistics. Physical Review E, 80(2), 026217.
"""
Base.@kwdef struct LMeasure{M, TM} <: AssociationMeasure
    K::Int = 2
    metric::M = Euclidean()
    tree_metric::TM = Euclidean()
    τx::Int = 1
    τy::Int = 1
    dx::Int = 2
    dy::Int = 2
    w::Int = 0
end

"""
    l_measure(measure::LMeasure, x::VectorOrDataset, y::VectorOrDataset)

Compute the [`LMeasure`](@ref) from source `x` to target `y`.
"""
function l_measure(measure::LMeasure, x::VectorOrDataset, y::VectorOrDataset)
    return estimate(measure, x, y)
end

function getrank(x, p)
    xmin, xmax = minimum(x), maximum
end

# Internal method for use with `independence`
function estimate(measure::LMeasure, x::AbstractDataset, y::AbstractDataset)
    (; K, metric, tree_metric, τx, τy, dx, dy, w) = measure

    # Match length of datasets by excluding end points.
    lx = length(x); ly = length(y)
    lx > ly ? X = x[1:ly, :] : X = x
    ly > lx ? Y = y[1:lx, :] : Y = y
    N = length(X)

    # Search for the K nearest neighbors of each points in both X and Y
    treeX = searchstructure(KDTree, X, tree_metric)
    treeY = searchstructure(KDTree, Y, tree_metric)
    neighborhoodtype, theiler = NeighborNumber(N), Theiler(0)
    idxs_X, dists_X = bulksearch(treeX, X, neighborhoodtype, theiler)
    idxs_Y = bulkisearch(treeY, Y, neighborhoodtype, theiler)
    Gᵢx = N / 2
    Gᵢᵏx = (K + 1) / 2
    Gᵢᵏxy = zeros(N)

    # Pre-allocate vectors to hold distances and ranks in inner loop
    dists_Gᵢᵏxy = zeros(K)
    ranks_Gᵢᵏxy = zeros(Int, K)
    for i in 1:N
        xᵢ = X[i]
        dxᵢ = dists_X[i]
        for j = 1:K
            wᵢⱼ = idxs_Y[i][j] # nearest neighbor indices in Y
            # TODO: get this distance from pre-computed distances? tricky indexing,
            # not sure if worth it.
            d_ᵢ_wᵢⱼ = evaluate(metric, xᵢ, X[wᵢⱼ])
            dists_Gᵢᵏxy[j] = d_ᵢ_wᵢⱼ
            # The distance is guaranteed to be within closed range of all other distances,
            # so no need to handle edge cases here.
            ranks_Gᵢᵏxy[j] = findfirst(d -> d ≈ d_ᵢ_wᵢⱼ, dxᵢ)
        end
        Gᵢᵏxy[i] = sum(ranks_Gᵢᵏxy) / K
    end

    return sum((Gᵢx .- Gᵢᵏxy) ./ (Gᵢx .- Gᵢᵏx)) / N
end
