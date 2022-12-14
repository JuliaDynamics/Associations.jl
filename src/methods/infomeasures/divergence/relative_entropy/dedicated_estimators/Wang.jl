export Wang

using Neighborhood: KDTree, Euclidean
using Neighborhood: bulksearch

"""
    Wang <: DivergenceEstimator
    Wang(k = 5, l = 5, w = 0)

The `Wang` relative entropy estimator (Wang et al., 2009[^Wang2009] computes the
relative entropy between two `d`-dimensional `Dataset`s using nearest neighbor
searches to estimate densities around points. `k` is the number of nearest neighbors
for the first dataset, and `l` is the number of neighbors for the second dataset.

`w` is the Theiler window, which controls how many temporal neighbors are excluded
during neighbor searches. `w = 0` means that only the point itself is excluded.

## Description

The `Wang` estimator computes the Shannon relative entropy, or KL-divergence, given by

```math
D(\\mathbb{P} || \\mathbb{Q}) = \\int_{\\mathbb{R}^d} p(x) \\log \\dfrac{p(x)}{q(x)} dx
```

where ``\\mathbb{P}`` and ``\\mathbb{Q}`` are continuous probability measures
with densities ``p(x)`` and ``q(x)`` (``x \\in \\mathcal{R}^D``),
with respect to the Lebesque measure ``\\mu``, and ``dx := \\mu(dx)``.

[^Wang2009]:
    Wang, Q., Kulkarni, S. R., & Verdú, S. (2009). Divergence estimation for
    multidimensional densities via k-Nearest-Neighbor distances. IEEE Transactions on
    Information Theory, 55(5), 2392-2405.
"""
Base.@kwdef struct Wang <: DivergenceEstimator
    k::Int = 5
    l::Int = 5
    w::Int = 0
end

function estimate(measure::RelativeEntropyRenyiDifferential, est::Wang,
        x::AbstractDataset{D},
        y::AbstractDataset{D}) where D
    q, base = measure.e.q, measure.e.base
    q == 1 || error("`divergence` not defined for `Wang` estimator for Rènyi with q = $(e.q)")
    (; k, l, w) = est
    n, m = length(x), length(y)
    tree_x = KDTree(x, Euclidean())
    tree_y = KDTree(y, Euclidean())
    # Eq. 17 in Wang allows for picking different k for each point in both spaces.
    # We here limit to fixed k and l
    ρs = last.(bulksearch(tree_x, x, NeighborNumber(k), Theiler(w))[2])
    νs = last.(bulksearch(tree_y, x, NeighborNumber(l))[2])

    hr = digamma(k) - digamma(l) + # correction
        (D / n) * sum(log.(νs ./ ρs)) +
        log(m / (n - 1))
    return hr / log(base, ℯ)
end
