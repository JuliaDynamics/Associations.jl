
export CMIEstimator
export condmutualinfo

"""
    ConditionalMutualInformation <: InformationMeasure
    CMI # alias

The supertype of all conditional mutual informations.
"""
abstract type ConditionalMutualInformation <: InformationMeasure end
const CMI = ConditionalMutualInformation

"""
The supertype of all conditional mutual information definitions.
"""
abstract type ConditionalMutualInformationDefinition <: Definition end
const CMIDefinition = ConditionalMutualInformationDefinition

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
    condmutualinfo(measure::CMI, est::ProbabilitiesEstimator, x, y)

    condmutualinfo(measure::CMI, est::CMIEstimator, x, y)

Estimate `measure`, a conditional mutual information of some kind, between `x` and `y` using the given
estimator.

- The first set of signatures is for discrete conditional mutual informations (meaning that they're
    computed directly from a set of probabilities estimated by a
    [`ProbabilitiesEstimator`](@ref)).
    For these methods to give meaningful results, you must ensure that the `est` gives
    the *the same* [`outcome_space`](@ref) for both `x` and `y`. For example, use a
    [`ValueHistogram`](@ref) with a [`FixedRectangularBinning`](@ref).
- The second set of signatures is for differential conditional mutual informations (meaning that they're
    computed from density estimates). For a full list of compatible definitions and
    estimators, see the online documentation.

Returns `div`, the divergence estimate, whose interpretation depends on the
combination of `definition` and `est`.

## Supported definitions

Cross-entropies appear in several forms in the literature. We currently support
the following measures:

- [`CMIShannon`](@ref). Discrete Shannon relative entropy.
"""
condmutualinfo(args...; kwargs...) = estimate(args...; kwargs...)

include("shannon/ConditionalMutualInformationShannon.jl")
include("renyi/ConditionalMutualInformationRenyi.jl")
include("estimators/estimators.jl")

# estimate(def::CMIDefinition, measure::CMI, est, x, y, z) =
#     estimate(def, measure, est, x, y, z)
# estimate(measure::CMI, est, x, y, z) =
#     estimate(measure, est, x, y, z)

# # Default to Shannon conditional mutual information
# estimate(est::ProbabilitiesEstimator, x, y, z) =
#     estimate(CMIShannon(), est, x, y, z)
