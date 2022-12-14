export DiscreteDivergence, DivergenceDefinition
"""
    DivergenceDefinition

The supertype for all divergence definitions.

## Interface

Subtypes are named `*Divergence_`, where `*` is the divergence type and  `_ = AuthorYEAR`.
Concrete subtypes are:

- [`TsallisDivergenceTsallis1998`](@ref)
- [`RenyiDivergenceRenyi1961`](@ref)

"""
abstract type DivergenceDefinition end

"""
    DiscreteDivergence <: DivergenceEstimator
    DiscreteDivergence(est::ProbabilitiesEstimators, definition::DivergenceDefinition)

`DiscreteDivergence` is a generic probability-based plug-in estimator for discrete
Shannon/Rényi/Tsallis relative entropy, and other divergences.

It is just a wrapper around a [`ProbabilitiesEstimator`](@ref), which controls how
probabilities are estimated from data, and a [`DivergenceDefinition`](@ref), which controls
which divergence formula these probabilities are plugged into.

## Compatible divergence definitions

- [`TsallisDivergenceTsallis1998`](@ref)
- [`RenyiDivergenceRenyi1961`](@ref)

## Compatible estimators

!!! warn "Common outcome space"
    Input variables must share the outcome space ``\\Omega`` for this definition
    to be valid. Not all [`ProbabilitiesEstimator`](@ref)s guarantee the same
    [`outcome_space`](@ref) when given different datasets. It is up to you to
    select a probabilities estimator that is meaningful in this context.

    Examples of meaningful estimators:
    - [`ValueHistogram`](@ref) with [`FixedRectangularBinning`](@ref), where the binning
        is independent on the input data.
"""
struct DiscreteDivergence{P<:ProbabilitiesEstimator, D <: DivergenceDefinition} <: DivergenceEstimator
    est::P
    definition::D
end

include("renyi.jl")
include("tsallis.jl")
