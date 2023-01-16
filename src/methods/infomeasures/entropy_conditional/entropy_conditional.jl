export entropy_conditional
export ConditionalEntropy
export ConditionalEntropyDefinition

"""
The supertype for all conditional entropies.
"""
abstract type ConditionalEntropy <: InformationMeasure end

"""
The supertype for all conditional entropy definitions.
"""
abstract type ConditionalEntropyDefinition <: Definition end

# Measures
include("CEShannon.jl")
include("CETsallisFuruichi.jl")
include("CETsallisAbe.jl")

entropy_conditional(measure::ConditionalEntropy, args...; kwargs...) =
    estimate(measure, args...; kwargs...)


"""
    entropy_conditional(measure::ConditionalEntropy, c::ContingencyMatrix{T, 2}) where T

Estimate the discrete version of the given [`ConditionalEntropy`](@ref) `measure` from
its direct (sum) definition, using the probabilities from a pre-computed
[`ContingencyMatrix`](@ref), constructed from two input variables `x` and `y`.

The convention is to compute the entropy of the variable in the *first* column of `c`
conditioned on the variable in the *second* column of `c`. To do the opposite, call this
function with a new contingency matrix where the order of the variables is reversed.

If `measure` is not given, then the default is `CEShannon()`.
"""
function entropy_conditional(measure::ConditionalEntropy, c::ContingencyMatrix)
    return estimate(measure, c)
end

"""
    entropy_conditional([measure::ConditionalEntropy], est::ProbabilitiesEstimator, x, y)

Estimate the conditional entropy `measure` between `x` and `y` by the difference of
the joint entropy and the marginal entropy of `y`, without any bias correction, using
the provided [`ProbabilitiesEstimator`](@ref) `est`.

Joint and marginal probabilities are computed by jointly discretizing `x` and `y` using
the approach given by `est`, and obtaining the marginal distribution from the joint
distribution.

This only works for estimators that have an implementation for
[`marginal_encodings`](@ref). See the
[online documentation](@ref probabilities_estimators_ce) for a list of
compatible measures.

If `measure` is not given, then the default is `CEShannon()`.
"""
function entropy_conditional(measure::ConditionalEntropy, est::ProbabilitiesEstimator, x, y)
    return estimate(measure, est, x, y)
end

"""
    entropy_conditional([measure::ConditionalEntropy], est::DifferentialEntropyEstimator, x, y)

Estimate the mutual information `measure` between `x` and `y` by a sum of three
entropy terms, without any bias correction, using any [`DifferentialEntropyEstimator`](@ref)
compatible with multivariate data.

See the [online documentation](@ref diffentropy_estimators_ce) for a list of
compatible measures.

If `measure` is not given, then the default is `CEShannon()`.
"""
function entropy_conditional(measure::ConditionalEntropy, est::DifferentialEntropyEstimator, x, y)
    return estimate(measure, est, x, y)
end

# Generic 3H-formulation of mutual information.
function marginal_entropies_ce2h(measure::ConditionalEntropy, est::ProbabilitiesEstimator, x, y)
    e = measure.e.definition
    @show "here"
    X = Dataset(x)
    Y = Dataset(y)
    XY = Dataset(X, Y)
    HY = entropy(e, est, Y)
    HXY = entropy(e, est, XY)
    return HY, HXY
end

# Override some definitions, because the estimator behaviour need to be adjusted
# for multiple input variables.
const WellDefinedCEProbEsts{m, D} = Union{
    SymbolicPermutation{m},
    ValueHistogram{<:FixedRectangularBinning{D}},
    ValueHistogram{<:RectangularBinning{T}},
    Dispersion
} where {m, D, T}

function marginal_entropies_ce2h(measure::ConditionalEntropy,
        est::WellDefinedCEProbEsts{m, D}, x, y) where {m, D}
    eX, eY = marginal_encodings(est, x, y)
    eXY = Dataset(eX, eY)
    e = measure.e
    HY = entropy(e, CountOccurrences(), eY)
    HXY = entropy(e, CountOccurrences(), eXY)
    return HY, HXY
end

function marginal_entropies_ce2h(measure::ConditionalEntropy, est::DifferentialEntropyEstimator, x, y)
    e = measure.e.definition
    X = Dataset(x)
    Y = Dataset(y)
    XY = Dataset(X, Y)
    HY = entropy(est, Y)
    HXY = entropy(est, XY)
    return HY, HXY
end
