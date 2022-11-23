export Wang

using Neighborhood: KDTree, Euclidean
using Neighborhood: bulksearch

"""
    Wang <: RelativeEntropyEstimator
    Wang(k = 1, w = 1)

The `Wang` relative entropy estimator (Wang et al., 2009[^Wang2009] computes the
relative entropy between two `d`-dimensional `Dataset`s using `k`-th nearest neighbor
searches to estimate densities around points.

[^Wang2009]:
    Wang, Q., Kulkarni, S. R., & Verdú, S. (2009). Divergence estimation for
    multidimensional densities via k-Nearest-Neighbor distances. IEEE Transactions on
    Information Theory, 55(5), 2392-2405.
"""
Base.@kwdef struct Wang <: RelativeEntropyEstimator
    k::Int = 1
    w::Int = 0
end

function entropy_relative(e::Renyi, est::Wang,
        x::AbstractDataset{D},
        y::AbstractDataset{D}) where D
    e.q == 1 || error("`entropy_relative` not defined for `Wang` estimator for Renyi with q = $(e.q)")
    (; k, w) = est
    n, m = length(x), length(y)

    tree_x = KDTree(x, Euclidean())
    tree_y = KDTree(y, Euclidean())
    # Eq. 17 in Wang allows for picking different k for each point in both spaces.
    # We here limit to fixed k.
    ds_x = last.(bulksearch(tree_x, x, NeighborNumber(k), Theiler(w))[2])
    ds_xiny = last.(bulksearch(tree_y, x, NeighborNumber(k), Theiler(w))[2])
    hr = (D / n) * sum(last.(ds_x) ./ last.(ds_xiny)) +
        log(m / (n - 1))
    return hr / log(e.base, ℯ)
end
entropy_relative(est::Wang, args...; base = 2, kwargs...) =
    entropy_relative(Shannon(; base), est, args...; kwargs...)
