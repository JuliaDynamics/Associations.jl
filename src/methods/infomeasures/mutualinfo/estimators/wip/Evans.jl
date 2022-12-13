using SpecialFunctions: gamma, digamma
using Neighborhood: bulksearch
using Neighborhood: KDTree, Euclidean, Theiler, NeighborNumber
using StateSpaceSets: Dataset
using StateSpaceSets: dimension

export Evans

"""
    Evans <: MutualInformationEstimator
    Evans(k = 1, w = 0)

The `Evans` mutual information estimator is based on `k`-th nearest neighbors
"""
Base.@kwdef struct Evans{M} <: MutualInformationEstimator
    k::Int = 5
    w::Int = 0
    metric::M = Chebyshev()
end

function estimate(def::MIShannonDifferential, est::Evans, x::Vector_or_Dataset...)
    e = def.e
    @assert length(x) >= 2 ||
        error("Need at leats two input datasets to compute mutual information between them.")

    (; k, w, metric) = est

    joint = Dataset(x...)
    marginals = Dataset.(x)
    N = length(joint)
    D = dimension(joint)
    tree_joint = KDTree(joint, metric)
    # The ball radii of interest are just the `k`-th nearest neighbor distances in the
    # joint and marginal spaces.
    rs_joint = last.(bulksearch(tree_joint, joint, NeighborNumber(k), Theiler(w))[2])
    rs = [zeros(N) for m in eachindex(marginals)]
    for (m, xₘ) in enumerate(marginals)
        distance_to_kth_neighbors!(est, rs[m], xₘ)
    end

    mi = 0.0
    marginal_rs = Dataset(rs...) # so we can index all marginals at once
    for (i, rs) in enumerate(marginal_rs)
        vⱼ = pball_volume(est, rs_joint[i]; d = D)
        vprod = prod(pball_volume(est, r; d = D) for r in rs)
        mi += log(vⱼ / vprod)
    end

    I = -digamma(k + 1) + digamma(N) - (mi / N)
    #@show -digamma(k)
    #@show digamma(N)
    #@show -(mi / N)
    return I / log(e.base, ℯ)
end
mutualinfo(est::Evans, args...; base = 2, kwargs...) =
    mutualinfo(Shannon(; base), est, args...; kwargs...)

# For the marginal dataset `xₘ`, find the distance to the `k`-th nearest neighbor
# of each point `xₘ[i]`, where `i = 1, 2, ..., N = length(xₘ)`.
function distance_to_kth_neighbors!(est::Evans, rs, xₘ)
    (; k, w, metric) = est
    tree = KDTree(xₘ, metric)
    rs[:] = last.(bulksearch(tree, xₘ, NeighborNumber(k), Theiler(w))[2])
end

# TODO: we could also implement this for Chebyshev distance.
"""
Compute the volume of a `d`-dimensional ball with radius `r` respecting the
Lₚ-norm.
"""
function pball_volume(p, r::Real = 1.0; d::Int)
    # https://link.springer.com/article/10.1007/s00013-019-01394-7
    if p == Inf
        return r^d
    end
    return (2*r)^d * gamma(1 + 1/p)^d / gamma(1 + d/p)
end

function pball_volume(est::Evans{M}, r::Real; d::Int) where M
    if M <: Euclidean
        p = 2
    elseif M <: Chebyshev
        p = Inf
    end
    return pball_volume(p, r; d)
end
