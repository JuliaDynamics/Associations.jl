using Entropies: ProbabilitiesEstimator, Entropy, EntropyEstimator, Shannon
using DelayEmbeddings: Dataset
export ConditionalMutualInformationEstimator
export cmi

"""  The supertype for all conditional mutual information definitions. """
abstract type ConditionalMutualInformation end

"""  The supertype for all conditional mutual information estimators """
abstract type ConditionalMutualInformationEstimator end

"""
    cmi(definition::ConditionalMutualInformation, est, x, y, z) → I::Real

Estimate the generalized conditional mutual information (CMI), using the formula and
logarithm base specified by `definition`) between `x` and `y` given `z`, using the
provided estimator.

## Supported definitions

Shannon's CMI is not the only possible definition of CMI. In the literature, there exists
multiple generalized CMI-like quantities that have been defined with respect to other
[`Entropy`](@ref) types. We currently support:

### Discrete

- [`CMIShannon`](@ref). Discrete Shannon CMI.

Discrete mutual information is computed when using an [`ProbabilitiesEstimator`](@ref).

### Continuous

- [`DifferentialCMIShannon`](@ref). Differential Shannon CMI.

Continuous mutual information is computed when using an [`EntropyEstimator`](@ref), a
[`MutualInformationEstimator`](@ref), or a dedicated
[`ConditionalMutualInformationEstimator`](@ref).

## Returns

`I` is the CMI, whose interpretation depends on the combination of `definition` and `est`.
"""
cmi(def::ConditionalMutualInformation, est, x, y) =
    estimate(def, est, x, y)

# Base.@kwdef struct CMI{METHOD} <: InformationMeasure
#     method::METHOD = nothing # e.g. 2MI, 4
# end

# """
#     cmi([e::Entropy,] est::ConditionalMutualInformationEstimator, x, y, z)
#     cmi([e::Entropy,] est::MutualInformationEstimator, x, y, z)
#     cmi([e::Entropy,] est::EntropyEstimator, x, y, z)
#     cmi([e::Entropy,] est::ProbabilitiesEstimator, x, y, z)

# Compute the discrete or continuous (differential) conditional mutual information (CMI)
# between `x` and `y` given `z`, i.e. ``II(X; Y | Z)``, using the provided estimator `est`.

# - If `est` is a [`ConditionalMutualInformationEstimator`](@ref), then the differential
#     CMI is computed using some direct estimation procedure.
# - If `est` is a [`MutualInformationEstimator`](@ref), then CMI is estimated by a sum of two
#     (differential) mutual information terms:
#     ``CMI(X; Y | Z) = I(X; Y, Z) + I(X; Z)``.
# - If `est` is an [`EntropyEstimator`](@ref), then, CMI is estimated by a sum of
#     four *differential entropy* terms
#     ``CMI(X; Y | Z) = h(X, Y) + h(Y, Z) - h(X, Y, Z) - h(Z)``.
# - If `est` is a [`ProbabilitiesEstimator`](@ref), then the discrete CMI is computed by
#     decomposing CMI into four discrete entropy terms:
#     ``CMI(X; Y | Z) = H(X, Y) + H(Y, Z) - H(X, Y, Z) - H(Z)``.

# If the entropy type (first argument) is not specified, then `Shannon(; base = 2)` is used.

# ## Common usage

# Used with [`LocalPermutation`](@ref) for conditional independence testing.


# """
# function cmi end

# This constant exist solely to allow nice default values. Add any
# new estimator types that are not `MutualInformationEstimator`s to this type union
# const CMI_ESTIMATOR_TYPES = Union{ProbabilitiesEstimator, EntropyEstimator, MutualInformationEstimator}


# function estimate(infomeasure::CMI{Nothing}, e::Entropy, est::CMI_ESTIMATOR_TYPES, x, y, z)
#     error("Please provide a valid estimation method to CMI, e.g. `CMI(MI2())`")
# end
# function estimate(infomeasure::CMI{MI2}, e::Entropy, est, x, y, z)
#     mutualinfo(e, est, x, Dataset(y, z)) -
#         mutualinfo(e, est, x, z)
# end
# function estimate(infomeasure::CMI{H4}, e::Entropy, est, x, y, z)
#     entropy(e, est, Dataset(x, z)) +
#         entropy(e, est, Dataset(y, z)) -
#         entropy(e, est, Dataset(x, y, z)) -
#         entropy(e, est, Dataset(z))
# end

# # Default to Shannon-type CMI
# cmi(est, x, y, z; base = 2) = cmi(Shannon(; base), est, x, y, z)
# # Default to using MI2 decomposition, because that'll work for both
# # EntropyEstimator/ProbabilityEstimator AND MutualInformationEstimator.
# # If using the H4 decomposition, MutualInformationEstimator won't work.
# # If feeding a
# cmi(e::Entropy, est, x, y, z;
#     method::EstimationMethod = MI2()) =
#     estimate(CMI(method), e, est, x, y, z)
# cmi(e::Entropy, est::ConditionalMutualInformationEstimator, x, y, z) =
#     estimate(CMI(), e, est, x, y, z)

# include("estimators/estimators.jl")
