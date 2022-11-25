using Neighborhood: Chebyshev, KDTree, NeighborNumber, Theiler
using Neighborhood: bulksearch
using SpecialFunctions: digamma
using DelayEmbeddings: AbstractDataset, Dataset
using DelayEmbeddings: dimension

export KraskovStögbauerGrassberger1, KSG1

"""
    KSG1 <: MutualInformationEstimator
    KraskovStögbauerGrassberger1 <: MutualInformationEstimator
    KraskovStögbauerGrassberger1(; k::Int = 1, w = 0, metric_marginals = Chebyshev())

The `KraskovStögbauerGrassberger1` mutual information estimator (you can use `KSG1` for
short) is the ``I^{(1)}`` `k`-th nearest neighbor estimator from
Kraskov et al. (2004)[^Kraskov2004].

## Keyword arguments

- **`k::Int`**: The number of nearest neighbors to consider. Only information about the
    `k`-th nearest neighbor is actually used.
- **`metric_marginals`**: The distance metric for the marginals for the marginals can be
    any metric from `Distances.jl`. It defaults to `metric_marginals = Chebyshev()`, which
    is the same as in Kraskov et al. (2004).
- **`w::Int`**: The Theiler window, which determines if temporal neighbors are excluded
    during neighbor searches in the joint space. Defaults to `0`, meaning that only the
    point itself is excluded.

## Description

Let the joint dataset ``X := \\{\\bf{X}_1, \\bf{X_2}, \\ldots, \\bf{X}_m \\}`` be defined by the
concatenation of the marginal datasets ``\\{ \\bf{X}_k \\}_{k=1}^m``, where each ``\\bf{X}_k``
is potentially multivariate. Let ``\\bf{x}_1, \\bf{x}_2, \\ldots, \\bf{x}_N`` be the points
in the joint space ``X``.

[^Kraskov2004]:
    Kraskov, A., Stögbauer, H., & Grassberger, P. (2004). Estimating mutual information.
    Physical review E, 69(6), 066138.
"""
struct KraskovStögbauerGrassberger1{MJ, MM} <: MutualInformationEstimator
    k::Int
    w::Int
    metric_joint::MJ # always Chebyshev, otherwise estimator is not valid!
    metric_marginals::MM

    function KraskovStögbauerGrassberger1(;
            k::Int = 1,
            w::Int = 0,
            metric_marginals::MM = Chebyshev()) where MM
        metric_joint = Chebyshev()
        new{typeof(metric_joint), MM}(k, w, metric_joint, metric_marginals)
    end
end

function estimate(infomeasure::MI{Nothing}, e::Renyi, est::KraskovStögbauerGrassberger1, x::Vector_or_Dataset...)
    @assert length(x) >= 2 ||
        error("Need at leats two input datasets to compute mutual information between them.")
    e.q == 1 || throw(ArgumentError(
        "Renyi entropy with q = $(e.q) not implemented for $(typeof(est)) estimators"
    ))
    (; k, w, metric_joint, metric_marginals) = est
    joint = Dataset(x...)
    marginals = Dataset.(x)
    M = length(x)
    N = length(joint)

    # `ds[i]` := the distance to k-th nearest neighbor for the point `joint[i]`.
    # In the paper, for method 1, ϵᵢ = 2*ds[i].
    tree_joint = KDTree(joint, metric_joint)
    ds = last.(bulksearch(tree_joint, joint, NeighborNumber(k), Theiler(w))[2])

    # `marginals_nₖs[i][k]` contains the number of points within radius `ds[i]` of
    # the point `marginals[k][i]`.
    ns = [zeros(Int, N) for m in eachindex(marginals)]
    s = 0.0
    for (m, xₘ) in enumerate(marginals)
        marginal_inrangecount!(est, ns[m], xₘ, ds)
    end
    marginal_nₖs = Dataset(ns...)

    mi = digamma(k) +
        (M - 1) * digamma(N) -
        (1 / N) * sum(sum(digamma.(nₖ)) for nₖ in marginal_nₖs)
    return mi / log(e.base, ℯ)
end
const KSG1 = KraskovStögbauerGrassberger1

function marginal_inrangecount!(est::KraskovStögbauerGrassberger1, ns, xₘ, ds)
    @assert length(ns) == length(xₘ)
    tree = KDTree(xₘ, est.metric_marginals)
    @inbounds for i in eachindex(xₘ)
        # Usually, we'd subtract 1 because inrangecount includes the point itself, but
        # we do the +1 in ψ(n + 1) here, so we don't need the +1 in the final MI formula.
        ns[i] = inrangecount(tree, xₘ[i], ds[i])
    end
    return ns
end
