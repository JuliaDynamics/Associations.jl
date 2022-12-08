using Neighborhood: KDTree, Euclidean, NeighborNumber, Theiler
using SpecialFunctions: digamma

export PoczosSchneiderRE

"""
    PoczosSchneiderRE <: RelativeEntropyEstimator
    PoczosSchneiderRE(k = 1, w = 0)

`PoczosSchneiderRE` is a relative entropy estimator that can be used to compute Renyi or
Tsallis relative entropies (divergences) (Póczos & Schneider, 2011)[^Póczos2011].

`w` is the Theiler window, which controls how many temporal neighbors are excluded
during neighbor searches. `w = 0` means that only the point itself is excluded.

## Description

Let ``\\mathbb{P}`` and ``\\mathbb{Q}`` be continuous probability measures
with density functions ``p(x)`` and ``q(x)`` (``x \\in M_0 \\subseteq \\mathcal{R}^D``),
i.e. ``p: M_0 \\subseteq \\mathcal{R}^D \\mapsto \\mathbb{R}`` and
``q: M_0 \\subseteq \\mathcal{R}^D \\mapsto \\mathbb{R}``,
with respect to the Lebesque measure ``\\mu``, and ``dx := \\mu(dx)``.

Let ``q\\in \\mathbb{R}`` with ``q \\neq 1``. Assuming the following
integral exists, `PoczosSchneiderRE` estimates

```math
V_{q}(\\mathbb{P} || \\mathbb{Q}) = \\int_{M_0} p^{q}(x)q^{1 - q}(x) dx.
```

The estimate for ``V_{q}(\\mathbb{P} || \\mathbb{Q})`` is then plugged into one of
the formulas below.

### Renyi divergence

If called with `entropy_relative(Renyi(), PoczosSchneiderRE(), x, y)`, then the Rényi
divergence is returned:

```math
R_{q}(\\mathbb{P} || \\mathbb{Q}) =
\\dfrac{1}{q - 1}\\log V_{q}(\\mathbb{P} || \\mathbb{Q}).
```

### Tsallis divergence

If called with `entropy_relative(Tsallis(), PoczosSchneiderRE(), x, y)`, then Tsallis
divergence is returned:

```math
T_{q}(\\mathbb{P} || \\mathbb{Q}) =
\\dfrac{1}{q - 1}
\\left(
    V_{q}(\\mathbb{P} || \\mathbb{Q}) - 1
\\right)
```

[^Póczos2011]:
    Póczos, B., & Schneider, J. (2011, June). On the estimation of alpha-divergences. In
    Proceedings of the Fourteenth International Conference on Artificial Intelligence and
    Statistics (pp. 609-617). JMLR Workshop and Conference Proceedings.
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
