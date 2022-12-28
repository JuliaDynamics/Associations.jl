export CrossEntropyDefinition
export DiscreteCrossEntropy

"""
The supertype of all cross entropy definitions.
"""
abstract type CrossEntropyDefinition end

"""
    DiscreteCrossEntropy <: CrossDifferentialEntropyEstimator
    DiscreteCrossEntropy(est::ProbabilitiesEstimator, definition::CrossEntropyDefinition)

`DiscreteCrossEntropy` is a generic probability-based plug-in estimator for discrete
generalized cross entropy.

It is just a wrapper around a [`ProbabilitiesEstimator`](@ref), which controls how
probabilities are estimated from data, and a [`CrossEntropyDefinition`](@ref), which
controls which cross entropy formula these probabilities are plugged into.

## Compatible definitions

Multiple generalizations of cross entropy exists, and `definiton` specifies which
generalization is to be computed. Currently, we support:

- [`RenyiDifferentialCrossEntropyThierrin2022`](@ref).

## Description

The `DiscreteCrossEntropy estimator first uses `est` to compute probabilities, then plugs
these probabilities into the cross entropy formula specified by `definition`.
"""
struct DiscreteCrossEntropy{P<:ProbabilitiesEstimator, D <: CrossEntropyDefinition} <: CrossDifferentialEntropyEstimator
    est::P
    definition::D
end

# TODO: maybe there should be a ContinuousCrossEntropy estimator too? (and the same for the other measures)

include("renyi/renyi.jl")
