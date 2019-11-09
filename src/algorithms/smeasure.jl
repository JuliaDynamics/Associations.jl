import Distances: Metric, SqEuclidean, evaluate, Euclidean
import NearestNeighbors: KDTree
import DynamicalSystems: FixedMassNeighborhood, neighborhood

"""
    s_measure(x, y, m::Int, τ, K::Int; distance_metric::Metric = SqEuclidean())

S-measure test [1] for the directional dependence between data series.

## Algorithm

1. Create an `m-`dimensional embeddings of both `x` and `y`, resulting in 
    `N` different `m`-dimensional embedding points
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
function s_measure(x, y, m::Int, τ, K::Int, metric = SqEuclidean(),
        tree_metric = Euclidean())
    
    length(x) == length(x) ? nothing : error("`x` and `y` must be same length")
    T = eltype(x)
    X = embed(x, m, τ)
    Y = embed(y, m, τ)
    N = length(X)
    
    treeX = KDTree(X, tree_metric)
    treeY = KDTree(Y, tree_metric)
    
    # Pre-allocate vectors to hold indices and distances during loops
    idxs_x = Vector{Int}(undef, K)
    idxs_y = Vector{Int}(undef, K)
    dists_x = Vector{T}(undef, K)
    dists_x_cond_y = Vector{T}(undef, K)
    
    # Mean squared distances in X
    Rx = Vector{T}(undef, N)
    
    # Mean squared distances in X conditioned on Y
    Rx_cond_y = Vector{T}(undef, N)
    
    neighborhoodtype = FixedMassNeighborhood(K)
    
    @inbounds for i in 1:N
        pxᵢ, pyᵢ = X[i], Y[i]
        
        idxs_x[:] = neighborhood(pxᵢ, treeX, neighborhoodtype, 1)
        idxs_y[:] = neighborhood(pyᵢ, treeY, neighborhoodtype, 1)
        
        for j = 1:K
            dists_x[j] = evaluate(metric, pxᵢ, X[idxs_x[j]])
            dists_x_cond_y[j] = evaluate(metric, pxᵢ, X[idxs_y[j]])
        end
        
        Rx[i] = sum(dists_x) / K
        Rx_cond_y[i] = sum(dists_x_cond_y) / K
    end
    
    S = sum(Rx ./ Rx_cond_y) / N
    return S
end

export s_measure