export BulinskiDimitrov

"""
    BulinskiDimitrov <: RelativeEntropyEstimator
    BulinskiDimitrov(k = 1, w = 0)

The `BulinskiDimitrov` estimator computes the [`Shannon`](@ref) differential cross
entropy between two variables based on `k`-th nearest neighbor searches,
using Theiler window `w`.
"""
Base.@kwdef struct BulinskiDimitrov{M} <: CrossEntropyEstimator
    k::Int = 1
    w::Int = 0
    metric::M = Euclidean()
end

function entropy_cross(e::Renyi, est::BulinskiDimitrov,
        x::AbstractDataset{D},
        y::AbstractDataset{D}) where D
    e.q == 1 || error("`entropy_cross` not defined for `BulinskiDimitrov` estimator for Renyi with q = $(e.q)")
    n, m = length(x), length(y)

    (; k, w, metric) = est
    n, m = length(x), length(y)
    tree_y = KDTree(y, metric)
    # Eq. 17 in Wang allows for picking different k for each point in both spaces.
    # We here limit to fixed k and l
    νs = last.(bulksearch(tree_y, x, NeighborNumber(k), Theiler(w))[2])

    # the first term is 1/N*sum(digamma.(kᵢ)) with `i = 1, 2, ..., n` if allowed to
    # vary.
    hc = digamma(k) +
        log(Entropies.ball_volume(D)) +
        log(m) +
        (D / n) * sum(log.(νs))

    return hc / log(e.base, ℯ)
end

entropy_cross(est::BulinskiDimitrov, args...; base = 2, kwargs...) =
    entropy_cross(Shannon(; base), est, args...; kwargs...)
