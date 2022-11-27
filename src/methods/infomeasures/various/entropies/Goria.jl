using Entropies: EntropyEstimator
using Entropies: ball_volume
using StateSpaceSets: AbstractDataset, Dataset
using Neighborhood: Euclidean, KDTree, NeighborNumber, Theiler
using Neighborhood: bulksearch
using SpecialFunctions: digamma

export Goria
"""
    Goria <: EntropyEstimator
    Goria(k = 1, w = 0, metric = Euclidean())

The `Goria` estimator estimates [`Shannon`](@ref) [`entropy`](@ref) by locating the
`k`-th nearest neighbors of each of the ``N`` query points, then
estimating entropy as described below (Goria et al., 2005)[^Goria2005].

## Description

Let ``\\bf{x}_1, \\bf{x}_2, \\ldots, \\bf{x}_N`` be input data points in ``\\mathbb{R}^m``,
and let ``\\bf{n}_1, \\bf{n}_2, \\ldots, \\bf{n}_N`` be the respective distances to their
`k`-th nearest neighbors. Next, let the geometric mean of the distances be

```math
\\hat{\\rho}_k = \\left( \\prod_{i=1}^N \\right)^{\\dfrac{1}{N}}
```

Goria et al. (2005)'s estimate of Shannon entropy is then

```math
\\hat{H} = m\\hat{\\rho}_k + \\log(N - 1) - \\psi(k) + \\log c_1(m),
```

where ``c_1(m) = \\dfrac{2\\pi^\\dfrac{m}{2}}{m \\Gamma(m/2)}`` and ``\\psi``
is the digamma function.

[^Goria2005]:
    Goria, M. N., Leonenko, N. N., Mergel, V. V., & Novi Inverardi, P. L. (2005). A new
    class of random vector entropy estimators and its applications in testing statistical
    hypotheses. Journal of Nonparametric Statistics, 17(3), 277-297.
"""
Base.@kwdef struct Goria{M} <: EntropyEstimator
    k::Int = 1
    w::Int = 0
    metric::M = Euclidean()
end

function entropy(e::Renyi, est::Goria, x::AbstractDataset{D}) where D
    e.q == 1 || throw(ArgumentError(
        "Renyi entropy with q = $(e.q) not implemented for $(typeof(est)) estimator"
    ))
    (; k, w, metric) = est
    N = length(x)

    tree = KDTree(x, metric)
    ds = last.(bulksearch(tree, x, NeighborNumber(k), Theiler(w))[2])
    h = D * log(prod(ds .^ (1 / N))) +
          log(N - 1) +
          log(c1(D)) -
          digamma(k)

    return h / log(e.base, ℯ)
end
c1(D::Int) = (2π^(D/2)) / (D* gamma(D/2))
