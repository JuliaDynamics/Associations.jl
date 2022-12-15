
export ConditionalMutualInformationEstimator
export cmi

"""
    ConditionalMutualInformation <: InformationMeasure

The supertype of all conditional mutual informations.
"""
abstract type ConditionalMutualInformation <: InformationMeasure end

"""
The supertype of all conditional mutual information definitions.
"""
abstract type ConditionalMutualInformationDefinition <: Definition end

"""
    ConditionalMutualInformationEstimator <: InformationEstimator

The supertype of all conditional mutual information estimators.

## Subtypes

- [``](@ref).
"""
abstract type ConditionalMutualInformationEstimator end

include("measures/ConditionalMutualInformationShannon.jl")
include("measures/ConditionalMutualInformationRenyi.jl")

include("definitions/CMI2MIShannon.jl")
include("definitions/CMI4HShannon.jl")

include("estimators/estimators.jl")


cmi(def::ConditionalMutualInformationDefinition, measure::ConditionalMutualInformation, est,
    x, y, z) = estimate(def, measure, est, x, y, z)
cmi(measure::ConditionalMutualInformation, est, x, y, z) = estimate(measure, est, x, y, z)

"""
    cmi(measure::ConditionalMutualInformation, est::ProbabilitiesEstimator, x, y)

    cmi(measure::ConditionalMutualInformation, est::ConditionalMutualInformationEstimator, x, y)

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
function cmi end

cmi(measure::ConditionalMutualInformation, est, x, y) = estimate(measure, est, x, y)
cmi(def::ConditionalMutualInformationDefinition, measure::ConditionalMutualInformation, est, x, y) =
    estimate(def, measure, est, x, y)
