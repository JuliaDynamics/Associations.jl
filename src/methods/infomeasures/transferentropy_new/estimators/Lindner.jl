using SpecialFunctions: digamma
using Neighborhood: KDTree, NeighborNumber, WithinRange, Euclidean
using Neighborhood: bulksearch, isearch

export Lindner

"""
    Lindner <: TransferEntropyEstimator
    Lindner(k = 1, w = 0, base = 2)

The `Lindner` transfer entropy estimator (Lindner et al., 2011)[^Lindner2011], which is
also used in the Trentool MATLAB toolbox, and is based on nearest neighbor searches.

`w` is the Theiler window, which determines if temporal neighbors are excluded
during neighbor searches (defaults to `0`, meaning that only the point itself is excluded
when searching for neighbours).

## Description

For a given points in the joint embedding space `jᵢ`, this estimator first computes the
distance `dᵢ` from `jᵢ` to its `k`-th nearest neighbor. Then, for each point `mₖ[i]` in
the `k`-th marginal space, it counts the number of points within radius `dᵢ`.

The transfer entropy is then computed as

```math
TE(X \\to Y) =
\\psi(k) + \\dfrac{1}{N} \\sum_{i}^n
\\left[
    \\sum_{k=1}^3 \\left( \\psi(m_k[i] + 1) \\right)
\\right],
```

where the index `k` references the three marginal subspaces `T`, `T𝒯` and `ST` for which
neighbor searches are performed.

[Lindner2011]: Lindner, M., Vicente, R., Priesemann, V., & Wibral, M. (2011). TRENTOOL:
    A Matlab open source toolbox to analyse information flow in time series data with
    transfer entropy. BMC neuroscience, 12(1), 1-22.
"""
Base.@kwdef struct Lindner{B} <: DifferentialEntropyEstimator
    k::Int = 2 # number of neighbors in joint space.
    w::Int = 0
    base::B = 2

    function Lindner(k::Int, w::Int, base::B) where B
        k >= 2 || throw(DomainError("The number of neighbors k must be >= 2."))
        new{B}(k, w, base)
    end
end

function transferentropy(measure::TEShannon, est::Lindner, args...)
    (; k, w, base) = est
    joint, ST, TT⁺, T = h4_marginals(measure, args...)
    N = length(joint)
    W = Theiler(w)
    metric =  Chebyshev()
    tree_joint = KDTree(joint, metric)
    nns_joint, ds_joint = bulksearch(tree_joint, joint, NeighborNumber(k), W)
    # For each `xᵢ ∈ M`, where `M` is one of the marginal spaces, count the number of
    # points within distance `ds[i]` from the point. Then count, for each point in each
    # of the marginals, how many neighbors each `xᵢ` has given `ds[i]`.
    ds = last.(ds_joint) # only care about distance to the k-th neighbor
    tree_ST = KDTree(ST, metric)
    tree_TT⁺ = KDTree(TT⁺, metric)
    tree_T = KDTree(T, metric)
    nns_ST  = [isearch(tree_ST, pᵢ, WithinRange(ds[i])) for (i, pᵢ) in enumerate(ST)]
    nns_TT⁺ = [isearch(tree_TT⁺, pᵢ, WithinRange(ds[i])) for (i, pᵢ) in enumerate(TT⁺)]
    nns_T   = [isearch(tree_T, pᵢ, WithinRange(ds[i])) for (i, pᵢ) in enumerate(T)]

    n_ST = length.(nns_ST)
    n_TT⁺ = length.(nns_TT⁺)
    n_T = length.(nns_T)
    te = 0.0
    for i = 1:N
        te += digamma(n_T[i] + 1) - digamma(n_TT⁺[i] + 1) - digamma(n_ST[i])
    end
    te /= N
    te += digamma(k)

    # Convert to target base *after* digamma computations, because the digamma function
    # is a function of the natural log.
    return te / log(ℯ, base)
end
