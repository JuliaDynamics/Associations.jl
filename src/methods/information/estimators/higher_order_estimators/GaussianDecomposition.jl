using LinearAlgebra: Symmetric, logdet
using Statistics: cov
using StatsBase: tiedrank
using Distributions: Normal, quantile

export GaussianDecomposition

"""
    GaussianDecomposition <: MultivariateInformationMeasureEstimator
    GaussianDecomposition(definition = OInformation(); copula = false)

A parametric estimator of [`HigherOrderInformationMeasure`](@ref)s that models the
inputs as jointly Gaussian and computes every joint entropy in closed form.

## Compatible definitions

- [`TotalCorrelation`](@ref)
- [`DualTotalCorrelation`](@ref)
- [`OInformation`](@ref)
- [`SInformation`](@ref)

## Description

The Shannon entropy of a ``d``-dimensional Gaussian with covariance matrix ``\\Sigma`` is

```math
H = \\frac{1}{2} \\log \\left( (2 \\pi e)^d \\det \\Sigma \\right).
```

The covariance matrix of all input columns is estimated once, and the entropy of every
subset of inputs uses the corresponding block of it. Inputs can be vectors or
multi-column `StateSpaceSet`s.

If `copula = true`, each column is first replaced by its Gaussian copula transform: its
ranks are mapped to standard normal quantiles, ``\\Phi^{-1}(r / (N + 1))``, following
[Ince2017](@cite). The estimate is then invariant to strictly monotonic transformations
of each column, and the estimator measures dependence in the Gaussian copula of the data.
In both cases the estimator only sees linear (or, with the copula, monotonic)
dependence: variables that are dependent but uncorrelated, such as ``x`` and ``x^2`` for a
symmetric ``x``, look independent to it.

No bias correction is applied.
"""
struct GaussianDecomposition{M<:MultivariateInformationMeasure} <: MultivariateInformationMeasureEstimator{M}
    definition::M
    copula::Bool
end

function GaussianDecomposition(definition = OInformation(); copula::Bool = false)
    return GaussianDecomposition(definition, copula)
end

function association(est::GaussianDecomposition{<:HigherOrderInformationMeasure}, x...)
    return higher_order_association(est, x...)
end

function subset_entropies(est::GaussianDecomposition, x...)
    X = StateSpaceSet.(x)
    widths = dimension.(X)
    M = reduce(hcat, Matrix.(X))
    if est.copula
        M = gaussian_copula_transform(M)
    end
    Σ = cov(M)
    ends = cumsum(widths)
    columns = [(e - w + 1):e for (e, w) in zip(ends, widths)]
    base = est.definition.base

    function h(S)
        idxs = reduce(vcat, collect.(columns[S]))
        d = length(idxs)
        hnats = 0.5 * (d * log(2π * ℯ) + logdet(Symmetric(Σ[idxs, idxs])))
        return convert_logunit(hnats, ℯ, base)
    end
    return SubsetEntropies(h)
end

# Replace each column by standard normal quantiles of its ranks.
function gaussian_copula_transform(M::AbstractMatrix)
    N = size(M, 1)
    out = similar(M, Float64)
    for (j, col) in enumerate(eachcol(M))
        out[:, j] .= quantile.(Normal(), tiedrank(col) ./ (N + 1))
    end
    return out
end
