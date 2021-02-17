
using Reexport

@reexport module SMeasure
export s_measure

import Distances: Metric, SqEuclidean, evaluate, Euclidean
import NearestNeighbors: KDTree
import ChaosTools: FixedMassNeighborhood, neighborhood
using Neighborhood
using ChaosTools
using DelayEmbeddings


"""
    s_measure(x::AbstractDataset, y::AbstractDataset; 
        K::Int = 2, metric = SqEuclidean(), tree_metric = Euclidean())
    s_measure(x::AbstractVector, y::AbstractVector; 
        K::Int = 2, metric = SqEuclidean(), tree_metric = Euclidean(),
        mx::Int = 2, my::Int = 2, τx::Int = 1, τy::Int = 1)
    s_measure(x::AbstractDataset, y::AbstractVector; 
        K::Int = 2, metric = SqEuclidean(), tree_metric = Euclidean(),
        mx::Int = 2, τx::Int = 1)
    s_measure(x::AbstractDataset, y::AbstractVector; 
        K::Int = 2, metric = SqEuclidean(), tree_metric = Euclidean(),
        my::Int = 2, τy::Int = 1)

S-measure test [1] for the directional dependence between time series `x` and `y`, which 
may be univariate or multivariate.

## Algorithm

1. If `x` and `y` are `m`-dimensional datasets, then use `x` and `y` as is. If `x` and `y` 
    are time series, then create an `m-`dimensional embeddings of both `x` and `y`, 
    resulting in `N` different `m`-dimensional embedding points
    ``X = \\{x_1, x_2, \\ldots, x_N \\}`` and ``X = \\{y_1, y_2, \\ldots, y_N \\}``.
    `τ` controls the embedding lag.
2. Let ``r_{i,j}`` and ``s_{i,j}`` be the indices of the `K`-th nearest neighbors 
    of ``x_i`` and ``y_i``, respectively.
3. Compute the the mean squared Euclidean distance to the ``K`` nearest neighbors 
    for each ``x_i``, using the indices ``r_{i, j}``.

```math
R_i^{(k)}(x) = \\dfrac{1}{k} \\sum_{i=1}^{k}(x_i, x_{r_{i,j}})^2
```

- Compute the y-conditioned mean squared Euclidean distance to the ``K`` nearest 
    neighbors for each ``x_i``, now using the indices ``s_{i,j}``.

```math
R_i^{(k)}(x|y) = \\dfrac{1}{k} \\sum_{i=1}^{k}(x_i, x_{s_{i,j}})^2
```

- Define the following measure of independence, where ``0 \\leq S \\leq 1``, and 
    low values indicate independence and values close to one occur for 
    synchronized signals.

```math
S^{(k)}(x|y) = \\dfrac{1}{N} \\sum_{i=1}^{N} \\dfrac{R_i^{(k)}(x)}{R_i^{(k)}(x|y)}
```

## References

Quian Quiroga, R., Arnhold, J. & Grassberger, P. [2000]
“Learning driver-response relationships from synchronization patterns,” 
Phys. Rev. E61(5), 5142–5148.
"""
function s_measure(x::AbstractDataset{D, T}, y::AbstractDataset{D, T}; K::Int = 3, 
        metric = SqEuclidean(),
        tree_metric = Euclidean(),
        theiler_window::Int = 0 # only point itself excluded
        ) where {D, T}

    length(x) == length(y) ? nothing : error("`x` and `y` must be same length")
    N = length(x)

    treeX = searchstructure(KDTree, x, tree_metric)
    treeY = searchstructure(KDTree, y, tree_metric)

    # Pre-allocate vectors to hold indices and distances during loops
    dists_x = Vector{T}(undef, K)
    dists_x_cond_y = Vector{T}(undef, K)

    # Mean squared distances in X
    Rx = Vector{T}(undef, N)

    # Mean squared distances in X conditioned on Y
    Rx_cond_y = Vector{T}(undef, N)

    # one more neighbor, because we need to excluded the first (itself) afterwards
    neighborhoodtype = NeighborNumber(K + 1) 
        
    idxs_X = bulkisearch(treeX, x, neighborhoodtype)
    idxs_Y = bulkisearch(treeY, y, neighborhoodtype)
    
    @inbounds for i in 1:N
        pxᵢ, pyᵢ = x[i], y[i]
    
        for j = 2:K
            dists_x[j] = @views evaluate(metric, pxᵢ, x[idxs_X[i][j]])
            dists_x_cond_y[j] = @views evaluate(metric, pxᵢ, x[idxs_Y[i][j]])
        end
        
        Rx[i] = sum(dists_x) / K
        Rx_cond_y[i] = sum(dists_x_cond_y) / K
    end

    S = sum(Rx ./ Rx_cond_y) / N
    return S
end

function s_measure(x::AbstractVector{T}, y::AbstractVector{T}; K::Int = 3,
    mx::Int = 2, my::Int = 2, τx = 1, τy = 1, metric = SqEuclidean(),
    tree_metric = Euclidean()) where T

    X = embed(x, mx, τx)
    Y = embed(y, my, τy)

    s_measure(X, Y, K = K, metric = metric, tree_metric = tree_metric)
end

function s_measure(x::AbstractDataset, y::AbstractVector{T}; K::Int = 3,
    my::Int = 2, τy = 1, 
    metric = SqEuclidean(), tree_metric = Euclidean()) where T

    Y = embed(y, my, τy)

    s_measure(x, Y, K = K, metric = metric, tree_metric = tree_metric)
end

function s_measure(x::AbstractVector{T}, y::AbstractDataset; K::Int = 3,
    mx::Int = 2, τx = 1, metric = SqEuclidean(),
    tree_metric = Euclidean()) where T

    X = embed(X, mx, τx)

    s_measure(X, y, K = K, metric = metric, tree_metric = tree_metric)
end


end