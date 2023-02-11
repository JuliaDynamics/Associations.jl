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
- Use with [`l_measure`](@ref) to compute the raw m-measure statistic.

## Description

[^Chicharro2009]:
    Chicharro, D., & Andrzejak, R. G. (2009). Reliable detection of directional couplings
    using rank statistics. Physical Review E, 80(2), 026217.
"""
Base.@kwdef struct LMeasure{M, TM} <: AssociationMeasure
    K::Int = 2
    metric::M = SqEuclidean()
    tree_metric::TM = Euclidean()
    τx::Int = 1
    τy::Int = 1
    dx::Int = 2
    dy::Int = 2
    w::Int = 0
end


function l_measure(measure::LMeasure, x::VectorOrDataset, y::VectorOrDataset)
    return estimate(measure, x, y)
end

# Internal method for use with `independence`
function estimate(measure::LMeasure, x::AbstractDataset, y::AbstractDataset)
    error("`estimate` not implemented for LMeasure")
    # (; K, metric, tree_metric, τx, τy, dx, dy, w) = measure

    # # Match length of datasets by excluding end points.
    # lx = length(x); ly = length(y)
    # lx > ly ? X = x[1:ly, :] : X = x
    # ly > lx ? Y = y[1:lx, :] : Y = y
    # N = length(X)

    # T = eltype(1.0)
    #  # Pre-allocate vectors to hold indices and distances during loops
    # dists_x = zeros(T, K)
    # dists_x_cond_y = zeros(T, K)

    # # Mean squared distances in X, and
    # # mean squared distances in X conditioned on Y
    # Rx = zeros(T, N)
    # Rx_cond_y = zeros(T, N)

    # # Search for the K nearest neighbors of each points in both X and Y
    # treeX = searchstructure(KDTree, X, tree_metric)
    # treeY = searchstructure(KDTree, Y, tree_metric)
    # neighborhoodtype, theiler = NeighborNumber(K), Theiler(w)
    # idxs_X = bulkisearch(treeX, X, neighborhoodtype, theiler)
    # idxs_Y = bulkisearch(treeY, Y, neighborhoodtype, theiler)

    # for n in 1:N
    #     pxₙ = X[n]

    #     for j = 1:K
    #         rₙⱼ = idxs_X[n][j] # nearest neighbor indices in X
    #         sₙⱼ = idxs_Y[n][j] # nearest neighbor indices in Y
    #         dists_x[j] = evaluate(metric, pxₙ, X[rₙⱼ])
    #         dists_x_cond_y[j] = evaluate(metric, pxₙ, X[sₙⱼ])
    #     end

    #     Rx[n] = sum(dists_x) / K
    #     Rx_cond_y[n] = sum(dists_x_cond_y) / K
    # end

    # return sum(log.((Rx .- Rx_cond_y^K) ./ (Rx .- Rx^K))) / N
end
