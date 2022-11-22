using Entropies: ProbabilitiesEstimator, Entropy, EntropyEstimator
export transferentropy, transferentropy!
using DelayEmbeddings: AbstractDataset

"""
Abstract type for transfer entropy estimators.

See also: [`MutualInformationEstimator`](@ref), [`EntropyEstimator`](@ref),
[`ProbabilitiesEstimator`](@ref).
"""
abstract type TransferEntropyEstimator end

# TODO: move all parameters here.
struct TransferEntropy <: InformationMeasure end

include("utils.jl")

function from_marginals(measure::TransferEntropy, e::Entropy, est::TransferEntropyEstimator,
    args...; kwargs...)
    msg = "$(typeof(e)) $(typeof(measure)) not implemented for transfer entropy estimator $(typeof(est))"
    throw(ArgumentError(msg))
end

function from_marginals(measure::TransferEntropy, e::Entropy, est::MutualInformationEstimator,
    args...; kwargs...)
    msg = "$(typeof(e)) $(typeof(measure)) not implemented for mutual information estimator $(typeof(est))"
    throw(ArgumentError(msg))
end

# function from_marginals(measure::TransferEntropy, e::Entropy, est::ProbabilitiesEstimator,
#         args...; kwargs...)
#     combo = "$(typeof(e)) transfer entropy, with probabilities estimated using $(typeof(est))"
#     msg = "Marginal estimation of $typeof(e) $(typeof(measure)) is not implemented for $combo)"
#     throw(ArgumentError(msg))
# end

VALID_TE_ESTIMATOR_TYPES = Union{
    EntropyEstimator,
    ProbabilitiesEstimator,
    MutualInformationEstimator,
    TransferEntropyEstimator}

"""
    transferentropy([e::Entropy,] est::ProbabilitiesEstimator, s, t, [c];
    transferentropy([e::Entropy,] est::EntropyEstimator, s, t, [c])
    transferentropy([e::Entropy,] est::MutualInformationEstimator, s, t, [c])
    transferentropy([e::Entropy,] est::TransferEntropyEstimator, s, t, [c])

Estimate transfer entropy of type `e` from the source timeseries `s` to the target
timeseries `t`, conditioned on timeseries `c` (if given).

The first argument - the entropy type - is optional and defaults to `Shannon()`.

## Estimation methods

Transfer entropy can be estimated using one of four estimator types

- **[`EntropyEstimator`](@ref)**. Decomposes the transfer entropy into a sum of marginal
    entropies, then computes the entropy of type `e` for each marginal.
    Example: [`Kraskov`](@ref).
- **[`ProbabilitiesEstimator`](@ref)**. Decomposes the transfer entropy into a sum of marginal
    entropies, then explicitly computes a probability distribution for each marginal based
    on some chosen property of the dat. Marginal entropies of type `e` are then computed by
    giving these probabilities to [`entropy`](@ref). Example: [`NaiveKernel`](@ref).
- **[`MutualInformationEstimator`](@ref)**. Decomposes the transfer entropy into a sum of
    mutual information terms, which are estimated separately using [`mutualinfo`](@ref).
    Example: [`Kraskov2`](@ref).
- **[`TransferEntropyEstimator`](@ref)**. Estimate the transfer entropy directly using the
    provided estimator.

Together, these estimators cover most existing transfer entropy estimation methods.
A (non-exhaustive) overview of implemented literature methods are found in the online
documentation.
"""
function transferentropy end
transferentropy(est::VALID_TE_ESTIMATOR_TYPES, args...; base = 2, kwargs...) =
    transferentropy(Shannon(; base), est, args...; kwargs...)


function transferentropy(e::Entropy, est::Union{EntropyEstimator, ProbabilitiesEstimator},
        args...; kwargs...)
    emb = EmbeddingTE(; kwargs...)
    joint, ST, T𝒯, T = get_marginals(TransferEntropy(), args...; emb)
    return entropy(e, est, T𝒯) +
        entropy(e, est, ST) -
        entropy(e, est, T) -
        entropy(e, est, joint)
end

# function from_marginals(m::TransferEntropy, e::Entropy, est::MutualInformationEstimator,
#     args...; kwargs...)
#     return
# end

# Estimators
include("estimators/estimators.jl")

# Literature methods.
include("convenience/symbolic/symbolic.jl")
include("convenience/hilbert/hilbert.jl")
include("convenience/transfer_operator/transferoperator.jl")
include("convenience/spike/spike.jl") # Continuous-time

# automated approaches
include("bbnue/bbnue.jl")
