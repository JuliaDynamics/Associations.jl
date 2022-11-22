using Entropies: ProbabilitiesEstimator, Entropy, EntropyEstimator, Shannon
using DelayEmbeddings: Dataset
export ConditionalMutualInformation
export conditional_mutualinfo

"""  The supertype for all conditional mutual information estimators """
abstract type ConditionalMutualInformationEstimator end

"""
    conditional_mutualinfo([e::Entropy,] est::ProbabilitiesEstimator, x, y, z)

Estimate ``I(x; y | z)``, the conditional mutual information between `x` and `y` given
`z`, by a sum of marginal entropies of type `e`, using the provided
[`ProbabilitiesEstimator`](@ref).

If the entropy type is not specified, then `Shannon(; base = 2)` is used.
"""
function conditional_mutualinfo(e::Entropy, est::ProbabilitiesEstimator,
        x::Vector_or_Dataset,
        y::Vector_or_Dataset,
        z::Vector_or_Dataset)
    mutualinfo(e, est, x, Dataset(y, z)) - mutualinfo(e, est, x, z)
end
conditional_mutualinfo(est::ProbabilitiesEstimator, args...; kwargs...) =
    conditional_mutualinfo(Shannon(; base = 2), est, args...; kwargs...)
conditional_mutualinfo(e::Entropy, x::Vector_or_Dataset, y::Vector_or_Dataset,
    z::Vector_or_Dataset) =
    error("Estimator missing. Please provide a valid estimator as the second argument.")

"""
    conditional_mutualinfo(e::EntropyEstimator, x, y, z)

Estimate ``I(x; y | z)``, the  conditional mutual information between `x` and `y` given
`z`, by a sum of marginal entropies (whose type is dictated by `e`), using the
provided [`EntropyEstimator`](@ref) estimator.
"""
function conditional_mutualinfo(e::EntropyEstimator,
        x::Vector_or_Dataset, y::Vector_or_Dataset, z::Vector_or_Dataset)
    mutualinfo(e, x, Dataset(y, z)) - mutualinfo(e, x, z)
end

"""
    conditional_mutualinfo(e::MutualInformation, x, y, z)
    conditional_mutualinfo(e::ConditionalMutualInformation, x, y, z)

Estimate ``I(x; y | z)`` by either

- using a dedicated [`MutualInformation`](@ref) estimator to compute
    ``I(x; y | z) = I(X; Y, Z) - I(X; Z)``.
- using a dedicated [`ConditionalMutualInformation`](@ref) estimator.

Both these ways of computing conditional mutual information may have less bias compared
to direct estimation using [`ProbabilityEstimator`](@ref)s or [`EntropyEstimator`](@ref)
estimators.

See also [`MutualInformation`](@ref) estimators [`Kraskov1`](@ref), [`Kraskov2`](@ref),
and [`ConditionalMutualInformation`](@ref) estimators.
"""
function conditional_mutualinfo(est::MutualInformationEstimator,
        x::Vector_or_Dataset, y::Vector_or_Dataset, z::Vector_or_Dataset)
    return mutualinfo(est, x, Dataset(y, z)) - mutualinfo(est, x, z)
end

function conditional_mutualinfo(est::ConditionalMutualInformationEstimator,
        x::Vector_or_Dataset, y::Vector_or_Dataset, z::Vector_or_Dataset)
    return mutualinfo(est, x, Dataset(y, z)) - mutualinfo(est, x, z)
end

include("estimators/estimators.jl")
