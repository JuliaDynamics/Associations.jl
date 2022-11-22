using Entropies: ProbabilitiesEstimator, EntropyEstimator, Entropy, Shannon

export MutualInformationEstimator
export mutualinfo

""" The supertype of all dedicated mutual information estimators """
abstract type MutualInformationEstimator end

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
I(X_1; X_2; \\ldots, X_m ) =
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

function mutualinfo(e::Entropy, est::ProbabilitiesEstimator,
        x::Vector_or_Dataset,
        y::Vector_or_Dataset)
    X = entropy(e, est, Dataset(x))
    Y = entropy(e, est, Dataset(y))
    XY = entropy(e, est, Dataset(x, y))
    MI = X + Y - XY
end
mutualinfo(est::ProbabilitiesEstimator, x::Vector_or_Dataset, y::Vector_or_Dataset) =
    mutualinfo(Shannon(; base = 2), est, x, y)
mutualinfo(e::Entropy, x::Vector_or_Dataset, y::Vector_or_Dataset) =
    error("Estimator missing. Please provide a valid estimator as the second argument.")

# """
#     mutualinfo([e::Entropy,] est::EntropyEstimator, x, y, ...)

# Estimate ``I(\\bf{X})``, the  mutual information between the datasets `x, z, ...`, by a
# sum of marginal entropies (whose type is dictated by `e`), using the provided
# [`EntropyEstimator`](@ref) estimator.
# """
function mutualinfo(e::Entropy, est::EntropyEstimator, x::Vector_or_Dataset...)
    @assert length(x) >= 2 ||
        error("Need at leats two input datasets to compute mutual information between them.")
    # Identity from Cover & Thomas (2006)
    h_joint = entropy(e, est, Dataset(x...))
    return sum(entropy(e, est, Dataset(xₖ)) for xₖ in x) - h_joint
end

mutualinfo(e::Entropy, est::EntropyEstimator, x::AbstractDataset) =
    mutualinfo(e, est, columns(x)...)

mutualinfo(est::Union{EntropyEstimator, ProbabilitiesEstimator}, args...;
    base = 2, kwargs...) =
    mutualinfo(Shannon(; base), est, args...; kwargs...)

include("estimators/estimators.jl")
