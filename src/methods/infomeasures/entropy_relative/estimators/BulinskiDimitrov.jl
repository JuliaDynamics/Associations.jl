export BulinskiDimitrov

Base.@kwdef struct BulinskiDimitrov{M} <: RelativeEntropyEstimator
    k::Int = 1
    w::Int = 0
    metric::M = Euclidean()
end

function entropy_relative(e::Renyi, est::BulinskiDimitrov,
        x::AbstractDataset{D},
        y::AbstractDataset{D}) where D
    e.q == 1 || error("`entropy_relative` not defined for `BulinskiDimitrov` estimator for Renyi entropy with q = $(e.q)")

    (; k, w, metric) = est
    n, m = length(x), length(y)
    tree_x = KDTree(x, metric)
    tree_y = KDTree(y, metric)

    # The equations allow variable k for each point in both space. We here fix it to k.
    Rs = last.(bulksearch(tree_x, x, NeighborNumber(k), Theiler(w))[2])
    Vs = last.(bulksearch(tree_y, x, NeighborNumber(k))[2])
    # The first two terms cancel when k is not variable.
    hr = digamma(k) - digamma(k) +
        log(m / (n - 1)) + # one point is excluded in the x-space
        D / n * sum(log.(Rs ./ Vs))

    return hr / log(e.base, ℯ)
end
entropy_relative(est::BulinskiDimitrov, args...; base = 2, kwargs...) =
    entropy_relative(Shannon(; base), est, args...; kwargs...)

function entropy(e::Renyi, est::BulinskiDimitrov, x::AbstractDataset{D}) where D
    e.q == 1 || error("`entropy` not defined for `BulinskiDimitrov` estimator for Renyi entropy with q = $(e.q)")

    (; k, w, metric) = est
    n = length(x)
    tree_x = KDTree(x, metric)
    # The equations allow variable k. We here fix it to k.
    Rs = last.(bulksearch(tree_x, x, NeighborNumber(k), Theiler(w))[2])
    # The first two terms cancel when k is not variable.
    logVd = log(Entropies.ball_volume(D))
    n1 = n - 1
    ψₖ = digamma(k)

    h = sum(log.((Rs .* logVd .* n1) ./ exp(ψₖ)))

    return (-h) / log(e.base, ℯ)
end
entropy(est::BulinskiDimitrov, args...; base = 2, kwargs...) =
    entropy(Shannon(; base), est, args...; kwargs...)
