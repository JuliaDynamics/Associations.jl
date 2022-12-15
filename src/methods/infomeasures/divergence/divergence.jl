export DivergenceEstimator
export DivergenceDefinition
export divergence

"""
The supertype for all divergences.
"""
abstract type Divergence <: InformationMeasure end

"""
The supertype for all divergence definitions.
"""
abstract type DivergenceDefinition <: Definition end

divergence(def::DivergenceDefinition, measure::Divergence, est, x, y) =
    estimate(def, measure, est, x, y)
divergence(measure::Divergence, est, x, y) =
    estimate(measure, est, x, y)

"""
    divergence(definition::RelativeEntropyShannon, est::ProbabilitiesEstimator, x, y) → div::Real
    divergence(definition::RelativeEntropyTsallis, est::ProbabilitiesEstimator, x, y) → div::Real
    divergence(definition::RelativeEntropyRenyi, est::ProbabilitiesEstimator, x, y) → div::Real

    divergence(definition::RelativeEntropyShannonDifferential, est::DivergenceEstimator, x, y) → div::Real
    divergence(definition::RelativeEntropyTsallisDifferential, est::DivergenceEstimator, x, y) → div::Real
    divergence(definition::RelativeEntropyRenyiDifferential, est::DivergenceEstimator, x, y) → div::Real

Estimate the divergence specified by `definition` (which also specifies the base of the
logarithm, if relevant) between two datasets `x` and `y`, using the provided estimator.

The first set of signatures is for discrete divergences (meaning that they're computed
directly from a set of probabilities estimated by a [`ProbabilitiesEstimator`](@ref)).
For these methods to give meaningful results, you must ensure that the `est` gives
the *the same* [`outcome_space`](@ref) for both `x` and `y`. For example, use a
[`ValueHistogram`](@ref) with a [`FixedRectangularBinning`](@ref).
The second set of signatures is for differential
divergences (meaning that they're computed from density estimates). For a full list of
compatible definitions and estimators, see the online documentation.

Returns `div`, the divergence estimate, whose interpretation depends on the
combination of `definition` and `est`.

## Supported definitions

Divergence measures are abundant in the literature. Sometimes, different authors give
different definitions of divergence measures with the same name. We currently support
the following measures (some of which may be tweaked according to multiple definitions):

- [`RelativeEntropyShannon`](@ref). Discrete Shannon relative entropy.
- [`RelativeEntropyRenyi`](@ref). Discrete Rényi relative entropy.
- [`RelativeEntropyTsallis`](@ref). Discrete Tsallis relative entropy.
- [`RelativeEntropyShannonDifferential`](@ref). Differential Shannon relative entropy.
- [`RelativeEntropyRenyiDifferential`](@ref). Differential Rényi relative entropy.
- [`RelativeEntropyTsallisDifferential`](@ref). Differential Tsallis relative entropy.
"""
function divergence end

divergence(measure, est, x, y) = estimate(measure, est, x, y)
# """
#     divergence(::Renyi, p::Probabilities, q::Probabilities)

# Estimate the (discrete) relative entropy, or KL divergence, between the pre-computed
# probability distributions `p` and `q`, where `p[i]` and `q[i]` is the probability of the
# `i`-th outcome in some [outcome_space](@ref) ``\\omega{X}``, defined as

# ```math
# D_{KL}(X || Y) = \\sum_{x \\in \\mathcal{X}} P(x) \\log{\\dfrac{P(x)}{Q(x)}}
# ```

# For this definition to be meaningful `p` and `q` must have *the same* outcome space.

# See also: [`probabilities`](@ref).
# """
function divergence end

include("relative_entropy/relative_entropy.jl")
