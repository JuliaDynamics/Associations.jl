export estimate
export InformationMeasure
export InformationMeasureEstimator
export InformationMeasureDefinition # Actually, InformationMeasure and InformationMeasureDefinition could be identical

"""
The supertype for all estimators of information-based measures.
"""
abstract type InformationMeasureEstimator end

"""
A generic supertype for definitions of information measures (one measure may have
multiple definitions).
"""
abstract type Definition end
"""
    InformationMeasure <: CausalityMeasure

The supertype for all definitions of information-based measures.

## Why use definitions?

Several information measures, such as mutual information, come in several forms
depending on what type of generalized entropy they are defined with respect to.
For example, there are at least three forms of Rényi mutual informations.

In CausalityTools.jl, each unique variant of a measure is a subtype of `InformationMeasure`.
For example, [`MITsallisFuruichi`](@ref) gives the formula for Furuichi (2006)'s
Rényi-based mutual information.

## Implemented measures

### Mutual informations

- [`MIShannon`](@ref). Discrete Shannon mutual information.
- [`MITsallisFuruichi`](@ref). Discrete Tsallis mutual information, as defined by
    Furuichi (2006).
"""
abstract type InformationMeasure <: CausalityMeasure end

"""
    estimate(e::Entropy, est::InformationMeasureEstimator, input::Vector_or_Dataset...)

Given some `input` data, estimate some information measure using the given
[`InformationMeasureEstimator`](@ref), with respect to the generalized entropy `e`.
"""
function estimate(measure::InformationMeasure, args...; kwargs...) end

include("discrete.jl")
include("mutualinfo/mutualinfo.jl")


# Things that will be eventually moved to Entropies.jl
include("various/probabilities.jl")
include("various/entropies.jl")
# """
#     EstimationMethod

# Some information measures, like mutual information or conditional mutual information,
# may be computed in sevaral different ways. In CausalityTools, we call these
# *compound* measures. For example, mutual information may be computed as a sum of
# marginal entropy terms, or in terms of a KL-divergence.

# `EstimationMethod` is the supertype of all compound estimation methods, which
# in the case of information measures are methods of decomposing higher-level measures
# into lower-level ones. Currently, subtypes are

# - [`MI2`](@ref) (used to estimate conditional mutual information)
# - [`H4`](@ref) (used to estimate conditional mutual information)
# - [`H3`](@ref) (used to estimate mutual information)
# """
# abstract type EstimationMethod end
# struct MI2 <: EstimationMethod end
# struct H3 <: EstimationMethod end
# struct H4 <: EstimationMethod end

# struct KLDiv <: EstimationMethod end

# export MI2, H3, H4

# Things that will be eventually moved to Entropies.jl
# include("various/probabilities.jl")
# include("various/entropies.jl")

# # The below belong in this package.
# include("entropy_cross/entropy_cross.jl")
# include("divergence/divergence.jl")
# include("entropy_conditional/entropy_conditional.jl")
# include("mutualinfo/mutualinfo.jl")
# include("conditional_mutualinfo/conditional_mutualinfo.jl")
#include("transferentropy/transferentropy.jl")
#include("predictive_asymmetry/predictive_asymmetry.jl")

# TODO:
# Explain somewhere in the documentation that if `ProbabilitiesEstimators`
# are used, then the *discrete* version of whatever measure is computed.
# Otherwise, the *differential* entropy/mi/cmi/whatever is estimated.
# Some methods are *compound measures*, in the sense that they can be built
# from lower-level measures.
