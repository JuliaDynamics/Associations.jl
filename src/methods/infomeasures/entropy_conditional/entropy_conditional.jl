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
This estimation method works for both numerical and categorical data.
If `measure` is not given, then the default is `CEShannon()`.

The convention is to compute the entropy of the variable in the *first* column of `c`
conditioned on the variable in the *second* column of `c`. To do the opposite, call this
function with a new contingency matrix where the order of the variables is reversed.

## Compatible measures

|                             | [`ContingencyMatrix`](@ref) |
| --------------------------- | :-------------------------: |
| [`CEShannon`](@ref)         |             ✓              |
| [`CETsallisFuruichi`](@ref) |             ✓              |
| [`CETsallisAbe`](@ref)      |             ✓              |
"""
function entropy_conditional(measure::ConditionalEntropy, c::ContingencyMatrix)
    return estimate(measure, c)
end

"""
    entropy_conditional([measure::ConditionalEntropy], est::ProbabilitiesEstimator, x, y)

Estimate the entropy of `x` conditioned on `y`, using the discrete version of the given
conditional entropy (CE) `measure`. The CE is computed the difference of
the joint entropy and the marginal entropy of `y`, using
the [`ProbabilitiesEstimator`](@ref) `est`, which must compatible with multivariate data
(that is, have an implementation for [`marginal_encodings`](@ref)).
No bias correction is applied. If `measure` is not given, then the default is `CEShannon()`.

## Estimators

Joint and marginal probabilities are computed by jointly discretizing `x` and `y` using
the approach given by `est`, and obtaining the marginal distribution for `y` from the joint
distribution.

| Estimator                    | Principle           | [`CEShannon`](@ref) | [`CETsallisAbe`](@ref) | [`CETsallisFuruichi`](@ref) |
| ---------------------------- | ------------------- | :-----------------: | :--------------------: | :-------------------------: |
| [`CountOccurrences`](@ref)   | Frequencies         |         ✓          |           ✓           |              x              |
| [`ValueHistogram`](@ref)     | Binning (histogram) |         ✓          |           ✓           |              x              |
| [`SymbolicPermutation`](@ref) | Ordinal patterns    |         ✓          |           ✓           |              x              |
| [`Dispersion`](@ref)         | Dispersion patterns |         ✓          |           ✓           |              x              |
"""
function entropy_conditional(measure::ConditionalEntropy, est::ProbabilitiesEstimator, x, y)
    return estimate(measure, est, x, y)
end

"""
    entropy_conditional([measure::ConditionalEntropy], est::DifferentialEntropyEstimator, x, y)

Estimate the entropy of `x` conditioned on `y`, using the differential/continuous
version of the given conditional entropy (CE) `measure`.  The CE is computed the difference of
the joint entropy and the marginal entropy of `y`, using
the [`DifferentialEntropyEstimator`](@ref) `est`, which must be compatible with multivariate data.
No bias correction is applied.
If `measure` is not given, then the default is `CEShannon()`.

## Estimators

| Estimator                        | Principle         | [`CEShannon`](@ref) | [`CETsallisAbe`](@ref) | [`CETsallisFuruichi`](@ref) |
| -------------------------------- | ----------------- | :-----------------: | :--------------------: | :-------------------------: |
| [`Kraskov`](@ref)                | Nearest neighbors |         ✓          |           x           |              x              |
| [`Zhu`](@ref)                    | Nearest neighbors |         ✓          |           x           |              x              |
| [`ZhuSingh`](@ref)               | Nearest neighbors |         ✓          |           x           |              x              |
| [`Gao`](@ref)                    | Nearest neighbors |         ✓          |           x           |              x              |
| [`Goria`](@ref)                  | Nearest neighbors |         ✓          |           x           |              x              |
| [`Lord`](@ref)                   | Nearest neighbors |         ✓          |           x           |              x              |
| [`LeonenkoProzantoSavani`](@ref) | Nearest neighbors |         ✓          |           x           |              x              |
"""
function entropy_conditional(measure::ConditionalEntropy, est::DifferentialEntropyEstimator, x, y)
    return estimate(measure, est, x, y)
end

# Generic 3H-formulation of mutual information.
function marginal_entropies_ce2h(measure::ConditionalEntropy, est::ProbabilitiesEstimator, x, y)
    e = measure.e.definition
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
