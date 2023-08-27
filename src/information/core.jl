# The "information" API relates to all measures that are direct functionals of
# (multivariate) probability mass functions or probability densities. This is
# the same approach as taken in ComplexityMeasures.jl.

export MultivariateInformationMeasure
export MultivariateInformationMeasureEstimator

const DiscreteOrDifferentialInfoEstimator = Union{OutcomeSpace, DiscreteInfoEstimator, DifferentialInfoEstimator}

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
abstract type MultivariateInformationMeasure <: AssociationMeasure end

abstract type MultivariateInformationMeasureEstimator end

"""
    estimate(e::InformationMeasure, est::InformationMeasureEstimator, input::VectorOrStateSpaceSet...)

Given some `input` data, estimate some information measure using the given
[`InformationMeasureEstimator`](@ref), with respect to the generalized entropy `e`.
"""
function estimate(measure::MultivariateInformationMeasure, args...; kwargs...) end
