using Statistics: quantile

import ComplexityMeasures: entropy
using LinearAlgebra: det
using SpecialFunctions: digamma
using Distributions: Distribution, UnivariateDistribution, quantile, MvNormal, Normal

export ParametricCopula

"""
    ParametricCopula <: MutualInformationEstimator
    ParametricCopula(d = Normal())

A parametric copula-based mutual information estimator.

Robin et al. (2016) A statistical framework for neuroimaging data analysis based on
mutual information estimated via a gaussian copula.
"""
Base.@kwdef struct ParametricCopula{D} <: MutualInformationEstimator
    d::D = Normal()
    debias = true
end

function estimate(measure::MIShannon, est::ParametricCopula, x, y)
    X = StateSpaceSet(x)
    Y = StateSpaceSet(y)
    D = StateSpaceSet(X, Y)
    Tp = StateSpaceSet((copula_transform(c) for c in columns(D))...)
    -entropy(est.d, Tp; debias = est.debias)
end

function entropy(d::Normal, x::AbstractStateSpaceSet{D}; debias = true, base = 2) where D
    N = length(x)
    Σ = fastcov(x)
    h = 1 / (2 * log(2))
    if debias
        bias = D * log(2/(N-1)) - sum(map(i -> digamma((N - i) / 2), 1:D))
        h *= log((2*π*ℯ)^D * det(Σ)) - bias
    else
        h *= log((2*π*ℯ)^D * det(Σ))
    end
    return _convert_logunit(h, ℯ, e.base)
end



"""
    inv_cdf_transform(x::AbstractVector, d::Distribution) → tx::Vector

Map each value `xᵢ ∈ x` to the transformed value `t(xᵢ) ∈ [0, 1]` using
the inverse cumulative distribution function (CDF) (i.e.e quantile function)
of the distribution `d`.

This function is meant to be used marginal empirical copula transforms (which
are uniformly distributed). Since the inverse CDF is a strictly increasing function,
the marginal copulas are preserved by the transformation.
"""
function inv_cdf_transform(x::AbstractVector, d::Distribution)
    ex = empcdf(x)
    t = zeros(length(ex))
    for (i, eᵢ) in enumerate(ex)
        if eᵢ == 1.0
            t[i] = quantile(d, 1-eps())
        elseif eᵢ == 0.0
            t[i] = quantile(d, eps())
        else
            t[i] = quantile(d, eᵢ)
        end
    end

    return t
end

function entropy_debiased(d::Normal, x::AbstractStateSpaceSet{D}) where D
    # Strictly speaking, `d` should be a MvNormal, but it doesn't matter here,
    # because the marginal data have already been transformed to be normally distributed.
    # `d` is purely for dispatch.
    N = length(x)
    Σ = fastcov(x)
    h = 1 / (2 * log(2))
    h *= log(2*π*ℯ^D * det(Σ)) - D * log(2/(N-1)) - sum(map(i -> digamma((N - i) / 2), 1:D))
end
