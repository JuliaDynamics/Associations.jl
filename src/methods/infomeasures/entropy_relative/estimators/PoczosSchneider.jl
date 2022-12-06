using Neighborhood: KDTree, Euclidean, NeighborNumber, Theiler
using SpecialFunctions: digamma

export PoczosSchneiderRE

"""
    PoczosSchneiderRE <: RelativeEntropyEstimator

A relative entropy estimator that can be used to compute Renyi or Tsallis
relative entropies (divergences) (Póczos & Schneider)[^Póczos2011].

## Definitions

We here re-state the definitions from Póczos & Schneider (2011).
Let ``p, q`` be density function ``\\mathbb{R}^d \\supseteq M_0 \\to \\mathbb{R}``
and let ``\\alpha \\in \\mathbb{R}`` with ``\\alpha \\neq 1``, and assume the following integrals
exist.

### Renyi divergence

```math
R_{\\alpha}(p || q) =
\\dfrac{1}{\\alpha - 1}\\log \\int_{M_0} p^{\\alpha}(x)q^{1-\\alpha}(x) dx
```

When `q = 1`, ``R_{\\alpha}(p || q) = KL(p || q)``.

### Tsallis divergence

```math
T_{\\alpha}(p || q) =
\\dfrac{1}{\\alpha - 1}
\\left(
    \\int_{M_0} p^{\\alpha}(x)q^{1-\\alpha}(x) dx - 1
\\right)
```

## Definitions

The integrals in both the Renyi and Tsallis divergences as stated here are identical.
The `PoczosSchneiderRE` estimator boils down to estimating this integral, then plugging the
estimate into the formulas above.

[^Póczos2011]:
    Póczos, B., & Schneider, J. (2011, June). On the estimation of alpha-divergences. In Proceedings of the Fourteenth International Conference on Artificial Intelligence and Statistics (pp. 609-617). JMLR Workshop and Conference Proceedings.
"""
Base.@kwdef struct PoczosSchneiderRE <: RelativeEntropyEstimator
    k::Int = 1
    w::Int = 0
end

function entropy_relative(e::Tsallis, est::PoczosSchneiderRE,
        x::AbstractDataset{D}, y::AbstractDataset{D}) where D
    D̂q = estimate_D̂(e.q, est, x, y)
    D̂qᵀ = 1 / (e.q - 1) * log(D̂q)
    return  D̂qᵀ / log(e.base, ℯ)
end

function entropy_relative(e::Renyi, est::PoczosSchneiderRE,
        x::AbstractDataset{D}, y::AbstractDataset{D}) where D
    D̂q = estimate_D̂(e.q, est, x, y)
    D̂qᴿ = 1 / (e.q - 1) * (D̂q - 1)
    return  D̂qᴿ / log(e.base, ℯ)
end

function estimate_D̂(q::Real, est::PoczosSchneiderRE,
        x::AbstractDataset{D}, y::AbstractDataset{D}) where D
    @assert est.k >= 2
    # q here is their α
    (; k, w) = est
    # Note: x and y are not required to have the same size.
    N, M = length(x), length(y)

    # Find neighbors
    tree_x = KDTree(x, Euclidean())
    tree_y = KDTree(y, Euclidean())
    idxs_x, ds_x = bulksearch(tree_x, x, NeighborNumber(k), Theiler(w))
    idxs_xiny, ds_xiny = bulksearch(tree_y, x, NeighborNumber(k), Theiler(w))
    ρs_x = last.(ds_x) .^ D
    νs_xiny = last.(ds_xiny) .^ D

    # Multiplicative bias and density factors can be pre-computed.
    Bₖ = gamma(k)^2 / (gamma(k - q + 1) * gamma(k + q - 1))
    f = (N - 1) / M * Bₖ
    @show Bₖ, f
    D̂ᵀ = 0.0
    for i = 1:N
        D̂ᵀ += (f * ρs_x[i] / νs_xiny[i])^(1 - q) * Bₖ
    end
    D̂ᵀ /= N
    return D̂ᵀ
end
