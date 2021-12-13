
using Reexport

@reexport module LMeasure
export l_measure

using NearestNeighbors, Distances, Neighborhood, DelayEmbeddings


"""
L-measure [^Chicharro2009] for the directional dependence between 
univariate/multivariate time series `x` and `y`. 

Returns a scalar `s ∈ [0, 1]`, where `s = 0` indicates independence between `x` and `y`, 
and higher values indicate synchronization between `x` and `y`, with complete 
synchronization for `s = 1.0`.
```

[^Chicharro2009]: Chicharro, D., & Andrzejak, R. G. (2009). Reliable detection of directional couplings using rank statistics. Physical Review E, 80(2), 026217.
"""
function l_measure(x::AbstractDataset{D1, T}, y::AbstractDataset{D2, T}; K::Int = 3, 
        metric = SqEuclidean(),
        tree_metric = Euclidean(),
        theiler_window::Int = 0 # only point itself excluded
        ) where {D1, D2, T}
    
    # Match length of datasets by excluding end points.
    lx = length(x); ly = length(y)
    lx > ly ? X = x[1:ly, :] : X = x
    ly > lx ? Y = y[1:lx, :] : Y = y
    N = length(X)
    
    # Pairwise distances between all points
    
    dists = Distances.colwise(Euclidean(), x, y)
    x = rand(100)
    y = rand(100)

    X, Y, dists = l_measure(x, y)
    #sp = zeros(Int, length(X)) 
    #for (i, xᵢ) in enumerate(X)
        #sortperm!(sp, dists[:, i]) 
    #end

    # vᵢⱼ := time indices of the *k* spatially nearest neighbors of `xᵢ`. vᵢⱼ[1] is the time index of the closest neighbor. 
    # wᵢⱼ := time indices of the *k* spatioanlly nearest neighbors of `yᵢ`. wᵢⱼ[1] is the time index of the closest neighbor. 
end

function l_measure(x::AbstractVector{T}, y::AbstractVector{T}; K::Int = 3,
    dx::Int = 2, dy::Int = 2, τx = 1, τy = 1, metric = SqEuclidean(),
    tree_metric = Euclidean()) where T

    X = embed(x, dx, τx)
    Y = embed(y, dy, τy)
    lX, lY = length(X), length(Y)

    l_measure(
        # TODO: cut the last points of the shortest resulting embedding.
        lX > lY ? X[1:lY, :] : X, 
        lY > lX ? Y[1:lX, :] : Y, 
        K = K, metric = metric, tree_metric = tree_metric)
end

end
