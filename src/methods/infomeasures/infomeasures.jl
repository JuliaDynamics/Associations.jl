export estimate
export InformationMeasure
export InformationMeasureEstimator
export InformationMeasureDefinition # Actually, InformationMeasure and InformationMeasureDefinition could be identical

const ProbOrDiffEst = Union{ProbabilitiesEstimator, DifferentialEntropyEstimator}
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
    InformationMeasure <: AssociationMeasure

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

### Conditional mutual information (CMI)

- [`CMIRenyiSarbu`](@ref). Discrete Rényi CMI.
"""
abstract type InformationMeasure <: AssociationMeasure end

"""
    estimate(e::EntropyDefinition, est::InformationMeasureEstimator, input::VectorOrStateSpaceSet...)

Given some `input` data, estimate some information measure using the given
[`InformationMeasureEstimator`](@ref), with respect to the generalized entropy `e`.
"""
function estimate(measure::InformationMeasure, args...; kwargs...) end

# Contingency matrices and its computation based on various probabilites
# estimators
include("marginal_encodings.jl")

# Things that will be eventually moved to ComplexityMeasures.jl
include("various/probabilities.jl")
include("various/entropies.jl")

# Higher-level measures
include("entropy_conditional/entropy_conditional.jl")
include("entropy_joint.jl")
include("mutualinfo/mutualinfo.jl")
include("condmutualinfo/condmutualinfo.jl")
include("transferentropy/transferentropy.jl")
include("predictive_asymmetry/predictive_asymmetry.jl") # old (TE-based)
include("predictive_asymmetry/PA.jl") # new
include("pmi.jl")
