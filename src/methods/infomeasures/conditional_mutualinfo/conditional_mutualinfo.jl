using Entropies: ProbabilitiesEstimator, Entropy, EntropyEstimator, Shannon
using DelayEmbeddings: Dataset
export ConditionalMutualInformationEstimator
export CMI
export cmi

"""  The supertype for all conditional mutual information estimators """
abstract type ConditionalMutualInformationEstimator end
Base.@kwdef struct CMI{METHOD} <: InformationMeasure
    method::METHOD = nothing # e.g. 2MI, 4
end

"""
    cmi([e::Entropy,] est::ProbabilitiesEstimator, x, y, z)

Compute the discrete conditional mutual information (CMI) betwen `x` and `y`, given `z`,
i.e. ``II(x; y | z)``, using the provided [`ProbabilitiesEstimator`](@ref). This
decomposes the CMI into four separate entropy terms, for which which probability
distributions are obtained, and then for which entropies are computed and summed.

    cmi([e::Entropy,] est::ConditionalMutualInformationEstimator, x, y, z)
    cmi([e::Entropy,] est::MutualInformationEstimator, x, y, z)
    cmi([e::Entropy,] est::EntropyEstimator, x, y, z)

Estimate the differential/continuous CMI using the provided estimators.
If `est` is a [`ConditionalMutualInformationEstimator`](@ref), then the differential
CMI is computed using some direct estimation procedure. If `est` is an
[`EntropyEstimator`](@ref), then the CMI is estimated by a sum of four marginal
differential entropies. If `est` is a [`MutualInformationEstimator`](@ref), then
the CMI is estimated by a sum of two (differential) mutual information terms.
"""
function cmi end

# """
#     conditional_mutualinfo([e::Entropy,] est::ProbabilitiesEstimator, x, y, z)
#     conditional_mutualinfo([e::Entropy,] est::EntropyEstimator, x, y, z)

# Estimate ``I(x; y | z)``, the conditional mutual information between `x` and `y` given
# `z`, by a sum of marginal entropies of type `e`, using the provided
# [`ProbabilitiesEstimator`](@ref).

# If the entropy type is not specified, then `Shannon(; base = 2)` is used.
# """
# function conditional_mutualinfo(e::Entropy, est::ProbabilitiesEstimator,
#         x::Vector_or_Dataset,
#         y::Vector_or_Dataset,
#         z::Vector_or_Dataset)
#     mutualinfo(e, est, x, Dataset(y, z)) - mutualinfo(e, est, x, z)
# end
# conditional_mutualinfo(est::ProbabilitiesEstimator, args...; kwargs...) =
#     conditional_mutualinfo(Shannon(; base = 2), est, args...; kwargs...)
# conditional_mutualinfo(e::Entropy, x::Vector_or_Dataset, y::Vector_or_Dataset,
#     z::Vector_or_Dataset) =
#     error("Estimator missing. Please provide a valid estimator as the second argument.")

# """
#     conditional_mutualinfo(e::EntropyEstimator, x, y, z)

# Estimate ``I(x; y | z)``, the  conditional mutual information between `x` and `y` given
# `z`, by a sum of marginal entropies (whose type is dictated by `e`), using the
# provided [`EntropyEstimator`](@ref) estimator.
# """


# This constant exist solely to allow nice default values. Add any
# new estimator types that are not `MutualInformationEstimator`s to this type union
const CMI_ESTIMATOR_TYPES = Union{ProbabilitiesEstimator, EntropyEstimator, MutualInformationEstimator}


function estimate(infomeasure::CMI{Nothing}, e::Entropy, est::CMI_ESTIMATOR_TYPES, x, y, z)
    error("Please provide a valid estimation method to CMI, e.g. `CMI(MI2())`")
end
function estimate(infomeasure::CMI{MI2}, e::Entropy, est, x, y, z)
    mutualinfo(e, est, x, Dataset(y, z)) -
        mutualinfo(e, est, x, z)
end

function estimate(infomeasure::CMI{H4}, e::Entropy, est, x, y, z)
    entropy(e, est, Dataset(x, y)) +
        entropy(e, est, Dataset(x, z)) -
        entropy(e, est, Dataset(z)) -
        entropy(e, est, Dataset(x, y, z))
end

# Default to Shannon-type CMI and estimating using the MI2 method
cmi(est, x, y, z; base = 2) = cmi(Shannon(; base), est, x, y, z)
cmi(e::Entropy, est, x, y, z; method::EstimationMethod = MI2()) =
    estimate(CMI(method), e, est, x, y, z)
cmi(e::Entropy, est::ConditionalMutualInformationEstimator, x, y, z) =
    estimate(CMI(), e, est, x, y, z)

# """
#     conditional_mutualinfo(e::MutualInformation, x, y, z)
#     conditional_mutualinfo(e::ConditionalMutualInformation, x, y, z)

# Estimate ``I(x; y | z)`` by either

# - using a dedicated [`MutualInformation`](@ref) estimator to compute
#     ``I(x; y | z) = I(X; Y, Z) - I(X; Z)``.
# - using a dedicated [`ConditionalMutualInformation`](@ref) estimator.

# Both these ways of computing conditional mutual information may have less bias compared
# to direct estimation using [`ProbabilityEstimator`](@ref)s or [`EntropyEstimator`](@ref)
# estimators.

# See also [`MutualInformation`](@ref) estimators [`Kraskov1`](@ref), [`Kraskov2`](@ref),
# and [`ConditionalMutualInformation`](@ref) estimators.
# """
# function conditional_mutualinfo(est::MutualInformationEstimator,
#         x::Vector_or_Dataset, y::Vector_or_Dataset, z::Vector_or_Dataset)
#     return mutualinfo(est, x, Dataset(y, z)) - mutualinfo(est, x, z)
# end

# function conditional_mutualinfo(est::ConditionalMutualInformationEstimator,
#         x::Vector_or_Dataset, y::Vector_or_Dataset, z::Vector_or_Dataset)
#     return mutualinfo(est, x, Dataset(y, z)) - mutualinfo(est, x, z)
# end

include("estimators/estimators.jl")
