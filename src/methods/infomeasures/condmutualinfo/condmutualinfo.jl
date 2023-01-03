
export CMIEstimator
export condmutualinfo

"""
    ConditionalMutualInformation <: InformationMeasure
    CMI # alias

The supertype of all conditional mutual informations.
"""
abstract type ConditionalMutualInformation{E, D} <: InformationMeasure end
const CMI{E, D} = ConditionalMutualInformation{E, D}

"""
The supertype of all conditional mutual information definitions.
"""
abstract type ConditionalMutualInformationDefinition <: Definition end
const CMIDefinition = ConditionalMutualInformationDefinition

""" The supertype for all H4-type (four entropies) decomposition of CMI. """
abstract type CMIH4 <: CMIDefinition end

"""
The supertype for all MI2-type (two mutual informations) decompositions of CMI.

Subtypes must implement the `measure::MutualInformation` with field, with a default value.
"""
abstract type CMIMI2{M <: MutualInformation} <: CMIDefinition end

"""
    ConditionalMutualInformationEstimator <: InformationEstimator
    CMIEstimator # alias

The supertype of all conditional mutual information estimators.

## Subtypes

- [`FrenzelPompeVelmejkaPalus`](@ref).
- [`PoczosSchneiderCMI`](@ref).
- [`Rahimzamani`](@ref).
- [`MesnerShalisi`](@ref).
"""
abstract type ConditionalMutualInformationEstimator end
const CMIEstimator = ConditionalMutualInformationEstimator


"""
    condmutualinfo(measure::CMI, est::CMIEstimator, x, y, z) → cmi::Real
    condmutualinfo(measure::CMI, est::DifferentialEntropyEstimator, x, y, z) → cmi::Real
    condmutualinfo(measure::CMI, est::ProbabilitiesEstimator, x, y, z) → cmi::Real

Estimate a conditional mutual information (CMI) of some kind (specified by `measure`),
between `x` and `y`, given `z`, using the given estimator.

## Definition

CMIs appear in many forms in the scientific literature. We support the following CMIs:

- **[`CMIShannon`](@ref)**: Shannon CMI.
- **[`CMIRenyi`](@ref)**: Renyi CMI.
"""
condmutualinfo(args...; kwargs...) = estimate(args...; kwargs...)

# - The first set of signatures is for discrete conditional mutual informations (meaning that
#     they're computed directly from a set of probabilities estimated by a
#     [`ProbabilitiesEstimator`](@ref)).
#     For these methods to give meaningful results, you must ensure that the `est` gives
#     the *the same* [`outcome_space`](@ref) for both `x` and `y`. For example, use a
#     [`ValueHistogram`](@ref) with a [`FixedRectangularBinning`](@ref).
# - The second set of signatures is for differential conditional mutual informations (meaning
#     that they're computed from density estimates). For a full list of compatible
#     definitions and estimators, see the online documentation.

include("definitions/definitions.jl")
include("CMIShannon.jl")
include("CMIRenyi.jl")

include("common_dispatch.jl")

include("estimators/estimators.jl")

# Default to Shannon mutual information.
condmutualinfo(est::ProbOrDiffEst, x, y, z) = estimate(CMIShannon(), est, x, y, z)
condmutualinfo(est::MutualInformationEstimator, x, y, z) =
    estimate(CMIShannon(), est, x, y, z)
