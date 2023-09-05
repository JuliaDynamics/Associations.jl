export MultivariateInformationMeasure
export MultivariateInformationMeasureEstimator

"""
    MultivariateInformationMeasure <: AssociationMeasure

The supertype for all multivariate information-based measure definitions.

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
