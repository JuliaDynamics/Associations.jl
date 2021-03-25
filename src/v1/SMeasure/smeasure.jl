
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
            K::Int = 2, metric = SqEuclidean(), tree_metric = Euclidean()) → Float64
        s_measure(x::AbstractVector, y::AbstractVector; 
            K::Int = 2, metric = SqEuclidean(), tree_metric = Euclidean(),
            dx::Int = 2, my::Int = 2, τx::Int = 1, τy::Int = 1) → Float64
        s_measure(x::AbstractDataset, y::AbstractVector; 
            K::Int = 2, metric = SqEuclidean(), tree_metric = Euclidean(),
            dx::Int = 2, τx::Int = 1) → Float64
        s_measure(x::AbstractDataset, y::AbstractVector; 
            K::Int = 2, metric = SqEuclidean(), tree_metric = Euclidean(),
            dy::Int = 2, τy::Int = 1) → Float64

    S-measure test [1] for the directional dependence between univariate/multivariate time 
    series `x` and `y`.

    ## Algorithm

    The algorithm is slightly modified from [1] to allow univariate time series as input.
    In that case, these time series are embedded using the parameters `mx`/`τx` (dimension/lag 
    for time series `x`) or `my`/`τy` (dimension/lag for time series `y`). For custom embeddings,
    do the embedding yourself and input `Dataset`s as both `x` and `y`.

    1a. If `x` and `y` are `dx` and `dy`-dimensional [`Dataset`](@ref)s, respectively, 
        then use `x` and `y` as is. 
    1b. If `x` and `y` are scalar time series, then create `dx` and `dy` dimensional embeddings,
        respectively, of both `x` and `y`, resulting in `N` different `m`-dimensional embedding points
        ``X = \\{x_1, x_2, \\ldots, x_N \\}`` and ``Y = \\{y_1, y_2, \\ldots, y_N \\}``.
        `τx` and `τy` control the embedding lags for `x` and `y`. 
        Datasets will be length-matched by eliminating points at the end of the longest dataset.
    1c. If `x` is a scalar-valued vector and `y` is a [`Dataset`](@ref), or vice versa, 
        then create an embedding of the scalar time series using parameters `dx`/`τx` or `dy`/`τy`.
        Datasets will be length-matched by eliminating points at the end of the longest dataset.
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
    function s_measure(x::AbstractDataset{D1, T}, y::AbstractDataset{D2, T}; K::Int = 3, 
            metric = SqEuclidean(),
            tree_metric = Euclidean(),
            theiler_window::Int = 0 # only point itself excluded
            ) where {D1, D2, T}
        
        size(x, 1) == size(y, 1) ? nothing : error("`x` and `y` must be same length")
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

        # use one more neighbor, because we need to excluded the first 
        # (itself) afterwards (distance to itself is zero)
        neighborhoodtype = NeighborNumber(K + 1) 
            
        idxs_X = bulkisearch(treeX, x, neighborhoodtype)
        idxs_Y = bulkisearch(treeY, y, neighborhoodtype)
        
        @inbounds for i in 1:N
            @views pxᵢ, pyᵢ = x[i], y[i]
        
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
        dx::Int = 2, dy::Int = 2, τx = 1, τy = 1, metric = SqEuclidean(),
        tree_metric = Euclidean()) where T

        X = embed(x, dx, τx)
        Y = embed(y, dy, τy)
        lX, lY = length(X), length(Y)

        s_measure(
            # TODO: cut the last points of the shortest resulting embedding.
            lX > lY ? X[1:lY, :] : X, 
            lY > lX ? Y[1:lX, :] : Y, 
            K = K, metric = metric, tree_metric = tree_metric)
    end

    function s_measure(X::AbstractDataset{D}, y::AbstractVector{T}; K::Int = 3,
            dy::Int = 2, τy = 1, 
            metric = SqEuclidean(), tree_metric = Euclidean()) where {D, T}

        Y = embed(y, dy, τy)
        s_measure(X[1:length(Y), :], Y, K = K, metric = metric, tree_metric = tree_metric)
    end

    function s_measure(x::AbstractVector{T}, Y::AbstractDataset{D}; K::Int = 3,
            dx::Int = 2, τx = 1, metric = SqEuclidean(),
            tree_metric = Euclidean()) where {D, T}
        X = embed(x, dx, τx)

        s_measure(X, Y[1:length(X), :], K = K, metric = metric, tree_metric = tree_metric)
    end


end