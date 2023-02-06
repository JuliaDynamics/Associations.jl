using NearestNeighbors, Distances, Neighborhood, DelayEmbeddings
using Distances: SqEuclidean, Euclidean

export SMeasure

"""
    SMeasure < AssociationMeasure
    SMeasure(K::Int = 2, metric = SqEuclidean(), tree_metric = Euclidean(),
        dx::Int = 2, my::Int = 2, τx::Int = 1, τy::Int = 1))

`SMeasure` is a bivariate association measure from Grassberger et al. (1999)[^Grassberger1999]
and Quiroga et al. (2000) [^Quiroga2000] that measure directional dependence 
between two input (potentially multivariate) time series.

## Usage

- Use with [`independence`](@ref) to perform a formal hypothesis test for directional dependence.
- Use with [`s_measure`](@ref) to compute the raw s-measure statistic.

## Description

The steps of the algorithm are:

1. Let ``r_{i,j}`` and ``s_{i,j}`` be the indices of the `K`-th nearest neighbors 
    of ``x_i`` and ``y_i``, respectively.

2. Compute the the mean squared Euclidean distance to the ``K`` nearest neighbors 
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

## Input data

The algorithm is slightly modified from [^Grassberger1999] to allow univariate timeseries as input.

- If `x` and `y` are [`Dataset`](@ref)s then use `x` and `y` as is and ignore the parameters
    `dx`/`τx` and `dy`/`τy`.
- If `x` and `y` are scalar time series, then create `dx` and `dy` dimensional embeddings,
    respectively, of both `x` and `y`, resulting in `N` different `m`-dimensional embedding points
    ``X = \\{x_1, x_2, \\ldots, x_N \\}`` and ``Y = \\{y_1, y_2, \\ldots, y_N \\}``.
    `τx` and `τy` control the embedding lags for `x` and `y`. 
- If `x` is a scalar-valued vector and `y` is a [`Dataset`](@ref), or vice versa, 
    then create an embedding of the scalar timeseries using parameters `dx`/`τx` or `dy`/`τy`.

In all three cases, input datasets are length-matched by eliminating points at the end of 
the longest dataset (after the embedding step, if relevant) before analysis.

[^Quiroga2000]:
    Quian Quiroga, R., Arnhold, J. & Grassberger, P. [2000] “Learning driver-response relationships
    from synchronization patterns,” Phys. Rev. E61(5), 5142–5148.
[^Grassberger1999]:
    Arnhold, J., Grassberger, P., Lehnertz, K., & Elger, C. E. (1999). A robust method for detecting
    interdependences: application to intracranially recorded EEG. Physica D:
    Nonlinear Phenomena, 134(4), 419-430.
"""
Base.@kwdef struct SMeasure{M, TM} <: AssociationMeasure
    K::Int = 2
    metric::M = SqEuclidean()
    tree_metric::TM = Euclidean()
    τx::Int = 1
    τy::Int = 1
    dx::Int = 2
    dy::Int = 2
end

function s_measure(measure::SMeasure, x::VectorOrDataset, y::VectorOrDataset)
    return estimate(measure, x, y)
end

# Internal method for use with `independence`
function estimate(measure::SMeasure, x::AbstractDataset, y::AbstractDataset)
    (; K, metric, tree_metric, τx, τy, dx, dy) = measure

    # Match length of datasets by excluding end points.
    lx = length(x); ly = length(y)
    lx > ly ? X = x[1:ly, :] : X = x
    ly > lx ? Y = y[1:lx, :] : Y = y
    N = length(X)

    T = eltype(1.0)
     # Pre-allocate vectors to hold indices and distances during loops
    dists_x = zeros(T, K)
    dists_x_cond_y = zeros(T, K)

    # Mean squared distances in X, and 
    # mean squared distances in X conditioned on Y
    Rx = zeros(T, N)
    Rx_cond_y = zeros(T, N)

    # Search for the K nearest neighbors of each points in both X and Y
    treeX = searchstructure(KDTree, X, tree_metric)
    treeY = searchstructure(KDTree, Y, tree_metric)
    neighborhoodtype, theiler = NeighborNumber(K), Theiler(0)
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
        
        Rx[n] = sum(dists_x) / K
        Rx_cond_y[n] = sum(dists_x_cond_y) / K
    end
    
    return sum(Rx ./ Rx_cond_y) / N
end

function estimate(measure::SMeasure, x::AbstractVector{T}, y::AbstractVector{T}) where T
    
    (; K, metric, tree_metric, τx, τy, dx, dy) = measure

    X = embed(x, dx, τx)
    Y = embed(y, dy, τy)
    lX, lY = length(X), length(Y)

    # TODO: cut the last points of the shortest resulting embedding.
    x̂ = lX > lY ? X[1:lY, :] : X
    ŷ = lY > lX ? Y[1:lX, :] : Y
    return estimate(measure, x̂, ŷ)
end

function estimate(measure::SMeasure, x::AbstractDataset{D}, y::AbstractVector{T}) where {D, T}
    (; K, metric, tree_metric, τx, τy, dx, dy) = measure

    Y = embed(y, dy, τy)
    X = x[1:length(Y), :]
    return estimate(measure, X, Y)
end

function estimate(measure::SMeasure, x::AbstractVector{T}, y::AbstractDataset{D}) where {D, T}
    (; K, metric, tree_metric, τx, τy, dx, dy) = measure

    X = embed(x, dx, τx)
    Y = y[1:length(X), :]
    return estimate(measure, X, Y)
end
