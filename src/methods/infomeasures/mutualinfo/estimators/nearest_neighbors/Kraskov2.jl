using Neighborhood: Chebyshev, KDTree, Theiler, NeighborNumber
using Neighborhood: bulksearch
using SpecialFunctions: digamma
using DelayEmbeddings: Dataset, AbstractDataset

export KraskovStögbauerGrassberger2, KSG2
"""
    KSG2 <: MutualInformationEstimator
    KraskovStögbauerGrassberger2 <: MutualInformationEstimator
    KraskovStögbauerGrassberger2(; k::Int = 1, w = 0, metric_marginals = Chebyshev(), base = 2)

The `KraskovStögbauerGrassberger2` mutual information estimator (you can use `KSG2` for
short) is the ``I^{(2)}`` `k`-th nearest neighbor estimator from
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

The `KraskovStögbauerGrassberger2` estimator first locates, for each ``\\bf{x}_i \\in X``, the
point ``\\bf{n}_i \\in X``, the `k`-th nearest neighbor to ``\\bf{x}_i``, according to the
maximum norm (`Chebyshev` metric). Let ``\\epsilon_i`` be the
distance ``d(\\bf{x}_i, \\bf{n}_i)``.

Consider ``x_i^m \\in \\bf{X}_m``, the ``i``-th point in the marginal space ``\\bf{X}_m``. For each
``\\bf{x}_i^m``, we determine ``\\theta_i^m`` := the number of points ``\\bf{x}_k^m \\in \\bf{X}_m`` that
are a distance less than ``\\epsilon_i`` away from ``\\bf{x}_i^m``. That is, we use the
distance from a query point ``\\bf{x}_i \\in X`` (in the *joint* space) to count neighbors of
``x_i^m \\in \\bf{X}_m`` (in the marginal space).

Mutual information between the variables ``\\bf{X}_1, \\bf{X_2}, \\ldots, \\bf{X}_m`` is
then estimated as

```math
\\hat{I}_{KSG2}(\\bf{X}) =
    \\psi{(k)} -
    \\dfrac{m - 1}{k} +
    (m - 1)\\psi{(N)} -
    \\dfrac{1}{N} \\sum_{i = 1}^N \\sum_{j = 1}^m \\psi{(\\theta_i^j + 1)}
```

[^Kraskov2004]:
    Kraskov, A., Stögbauer, H., & Grassberger, P. (2004). Estimating mutual information.
    Physical review E, 69(6), 066138.
"""
struct KraskovStögbauerGrassberger2{MJ, MM, B} <: MutualInformationEstimator
    k::Int
    w::Int
    metric_joint::MJ # always Chebyshev, otherwise estimator is not valid!
    metric_marginals::MM
    base::B

    function KraskovStögbauerGrassberger2(;
            k::Int = 1,
            w::Int = 0,
            metric_marginals::MM = Chebyshev(),
            base::B = 2) where {MJ, MM, B}
        metric_joint = Chebyshev()
        new{typeof(metric_joint), MM, B}(k, w, metric_joint, metric_marginals, base)
    end
end

const KSG2 = KraskovStögbauerGrassberger2

function mutualinfo(e::Renyi, est::KraskovStögbauerGrassberger2, x::Vector_or_Dataset...)
    @assert length(x) >= 2 ||
        error("Need at leats two input datasets to compute mutual information between them.")
    e.q == 1 || throw(ArgumentError(
        "Renyi entropy with q = $(e.q) not implemented for $(typeof(est)) estimators"
    ))
    (; k, metric_joint, metric_marginals, base) = est
    joint = Dataset(x...)
    marginals = Dataset.(x)
    M = length(x)
    N = length(joint)
    D = dimension(joint)

    # `ds[i]` := the distance to k-th nearest neighbor for the point `joint[i]`.
    tree_joint = KDTree(joint, metric_joint)
    knn_idxs = last.(bulksearch(tree_joint, joint, NeighborNumber(k), Theiler(0))[1])
    ns = [zeros(Int, N) for _ in eachindex(marginals)]
    s = 0.0
    for (k, xᵢ) in enumerate(marginals)
        marginal_inrangecount!(est, ns[k], xᵢ, knn_idxs)
    end
    marginal_nₖs = Dataset(ns...)

    h = digamma(k) -
        (M - 1) / k +
        (M - 1) * digamma(N) -
        (1 / N) * sum(sum(digamma.(nₖ)) for nₖ in marginal_nₖs)
    return h / log(base, ℯ)
end
mutualinfo(est::KraskovStögbauerGrassberger2, args...) = mutualinfo(Shannon(), est, args...)

function marginal_inrangecount!(est::KraskovStögbauerGrassberger2, ns::Vector{Int},
        xₘ::AbstractDataset, knn_idxs)
    (; k, w, metric_joint, metric_marginals, base) = est
    @assert length(ns) == length(xₘ)

    tree = KDTree(xₘ, metric_marginals)
    @inbounds for (i, (xᵢᵐ, joint_nnᵢ)) in enumerate(zip(xₘ, knn_idxs))
        # ϵᵢᵐ := distance from marginal point xᵢᵐ to the joint space query point `joint[i]`,
        # considering only the marginal axes for `xₘ`.
        ϵᵢᵐ = evaluate(Chebyshev(), xᵢᵐ, xₘ[joint_nnᵢ])

        # Usually, we'd subtract 1 because inrangecount includes the point itself, but
        # we do the +1 in ψ(n + 1) here, so we don't need the +1 in the final MI formula.
        ns[i] = inrangecount(tree, xᵢᵐ, ϵᵢᵐ)
    end
    return ns
end
