using Entropies: ProbabilitiesEstimator, EntropyEstimator, Entropy, Shannon

export MutualInformationEstimator
export MI
export mi

""" The supertype of all dedicated mutual information estimators """
abstract type MutualInformationEstimator end

Base.@kwdef struct MI{METHOD} <: InformationMeasure
    method::METHOD = nothing # e.g. 3H
end

"""
    mutualinfo([e::Entropy,] est::ProbabilitiesEstimator, X₁, X₂, ...)
    mutualinfo([e::Entropy,] est::EntropyEstimator, X₁, X₂, ...)
    mutualinfo([e::Entropy,] est::MutualInformationEstimator, X₁, X₂, ...)

Estimate ``I(\\bf{X})``, the  mutual information between the datasets
``\\bf{X} = \\{\\bf{X}_1,\\bf{X}_2, \\ldots, \\bf{X}_m \\}``, using the provided
estimator `est`, where `b := est.base` specifies the logarithm.

The first argument, the entropy type `e`, is optional and defaults to `Shannon()`.
Mutual information may also be defined for other entropy types, such as [`Renyi`](@ref)
or [`Tsallis`](@ref).

## Description

Mutual information ``I`` between ``X`` and ``Y`` is defined as

```math
I(X_1; X_2; \\ldots; X_m ) =
\\sum_{x_1, x_2, \\ldots, x_m}
p(x_1, x_2, \\ldots, x_m)
\\log_b \\left( \\dfrac{p(x_1, x_2, \\ldots, x_m)}{p(x_1)p(x_2)\\cdots p(x_m)} \\right)
```

If `est` is a [`MutualInformationEstimator`](@ref), then the mutual information is computed
using some specialized algorithm. If `est` is a [`ProbabilitiesEstimator`](@ref) or
[`EntropyEstimator`](@ref), then ``I(\\bf{X})`` is estimated as

```math
\\hat{I}(\\bf{X}) = - H_e(\\bf{X}) + \\sum_{i=1}^k H_e(X_i)
```

where ``H_e(\\cdot)`` is the entropy of type `e`.

[^Cover2006]: Cover, T. M., Thomas, J. A. (2006). Elements of Information Theory 2nd
    Edition (Wiley Series in Telecommunications and Signal Processing). Wiley-Interscience.
    ISBN: 0471241954
"""
function mutualinfo end


# This constant exist solely to allow nice default values. Add any
# new estimator types that are not `MutualInformationEstimator`s to this type union
const MI_ESTIMATOR_TYPES = Union{ProbabilitiesEstimator, EntropyEstimator}

# If `x` is variable, then H3 is treated as a N-component estimate
# We just call it H3, because bivariate MI it is the most common use case.
function estimate(infomeasure::MI{H3}, e::Entropy, est::MI_ESTIMATOR_TYPES,
        x::Vector_or_Dataset...)
    @assert length(x) >= 2 ||
        error("Need at leats two input datasets to compute mutual information between them.")
    # Identity from Cover & Thomas (2006)
    h_joint = entropy(e, est, Dataset(x...))
    return sum(entropy(e, est, Dataset(xₖ)) for xₖ in x) - h_joint
end

# Informative error messages.
function estimate(infomeasure::MI{Nothing}, e::Entropy, est::MI_ESTIMATOR_TYPES, args...)
    error("""Please provide a valid estimation method to MI, e.g.
    `estimate(MI(H3()), Shannon(), KSG1(), x, y, ...)`""")
end
function estimate(infomeasure::MI{Nothing}, e::Entropy, est::MutualInformationEstimator, args...;
        kwargs...)
    error("`estimate` not implemented for `MI` with estimator $(typeof(est))")
end
function estimate(infomeasure::MI{Nothing}, est::MutualInformationEstimator, args...;
        kwargs...)
    error("""Entropy type `e` must be specified, e.g.
    `estimate(MI(), Shannon(), KSG1(), x, y)}`""")
end

# Default to Shannon-type MI and estimating using the H3 method
# If dedicated estimators have other defaults, override in `./estimators/relevant_file.jl`.
mutualinfo(est::MI_ESTIMATOR_TYPES, x...; base = 2) =
    mutualinfo(Shannon(; base), est, x...)
mutualinfo(e::Entropy, est::MI_ESTIMATOR_TYPES, x...; method::EstimationMethod = H3()) =
    estimate(MI(method), e, est, x...)
mutualinfo(e::Entropy, est::MutualInformationEstimator, x...) =
    estimate(MI(), e, est, x...)

include("estimators/estimators.jl")
