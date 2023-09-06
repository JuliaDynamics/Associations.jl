export MutualInformation
export MutualInformationEstimator
export MutualInformationDefinition
export mutualinfo

"""
    MutualInformation <: AssociationMeasure

The supertype of all mutual information measures. Concrete subtypes are

- [`MIShannon`](@ref)
- [`MITsallisFuruichi`](@ref)
- [`MITsallisMartin`](@ref)
- [`MIRenyiJizba`](@ref)
- [`MIRenyiSarbu`](@ref)
"""
abstract type MutualInformation{E} <: InformationMeasure end

"""
    MutualInformationEstimator

The supertype of all dedicated mutual information estimators.

[`MutualInformationEstimator`](@ref)s can be either mixed, discrete or a combination of
both. Each estimator uses a specialized technique to approximate relevant
densities/integrals and/or probabilities, and is typically tailored to a specific
type of [`MutualInformation`](@ref) (mostly [`MIShannon`](@ref)).
"""
abstract type MutualInformationEstimator <: InformationMeasureEstimator end

# There are many ways of defining mutual information. Moreover, the definitions
# differ for different types of base `EntropyDefinition`s. Therefore, we dispatch
# on subtypes of `MutualInformationDefinition`.
""" The supertype for mutual information definitions. """
abstract type MutualInformationDefinition <: Definition end

""" The supertype for all H3-type (three entropies) decomposition of mutual information. """
abstract type MIH3 <: MutualInformationDefinition end

#= """
    mutualinfo(measure::MutualInformation, est::MutualInformationEstimator, x, y)
    mutualinfo(measure::MutualInformation, est::DifferentialEntropyEstimator, x, y)
    mutualinfo(measure::MutualInformation, est::ProbabilitiesEstimator, x, y)
    mutualinfo(measure::MutualInformation, c::ContingencyMatrix)

Estimate the mutual information `measure` (either [`MIShannon`](@ref) or
[`MITsallis`](@ref), ) between `x` and `y` using the provided estimator `est`.
Alternatively, compute mutual information from a pre-computed [`ContingencyMatrix`](@ref).

Compatible measures/definitions and estimators are listed in the
[online documentation](@ref mutualinfo_overview).
""" =#
mutualinfo(args...; kwargs...) = estimate(args...; kwargs...)

include("MIShannon.jl")
include("MITsallisFuruichi.jl")
include("MITsallisMartin.jl")
include("MIRenyiSarbu.jl")
include("MIRenyiJizba.jl")

include("estimators/estimators.jl")

# Default to Shannon mutual information.

"""
    mutualinfo([measure::MutualInformation], m::ContingencyMatrix) → mi::Real

Estimate the mutual information between `x` and `y`, the variables corresponding to
the columns and rows of the 2-dimensional contingency matrix `m`, respectively.

Estimates the discrete version of the given [`MutualInformation`](@ref) `measure` from
its direct definition (double-sum), using the probabilities from a pre-computed
[`ContingencyMatrix`](@ref). If `measure` is not given, then the default
is `MIShannon()`.
"""
mutualinfo(c::ContingencyMatrix) = estimate(MIShannon(), c)

"""
    mutualinfo([measure::MutualInformation], est::ProbabilitiesEstimator, x, y) → mi::Real ∈ [0, a]

Estimate the mutual information between `x` and `y` using the discrete version of the
given `measure`, using the given [`ProbabilitiesEstimator`](@ref) `est` (which must accept
multivariate data and have an implementation for [`marginal_encodings`](@ref)).
See examples [here](@ref example_mi_ProbabilitiesEstimator).
If `measure` is not given, then the default is `MIShannon()`.

## Estimators

The mutual information is computed as sum of three entropy terms, without any bias correction.
The exception is when using [`Contingency`](@ref); then the mutual information
is computed using a [`ContingencyMatrix`](@ref).

Joint and marginal probabilities are computed by jointly discretizing `x` and `y` using
the approach given by `est` (using [`marginal_encodings`](@ref)), and obtaining marginal
distributions from the joint distribution.

| Estimator                 | Principle           | [`MIShannon`](@ref) | [`MITsallisFuruichi`](@ref) | [`MITsallisMartin`](@ref) | [`MIRenyiJizba`](@ref) | [`MIRenyiSarbu`](@ref) |
| ------------------------- | ------------------- | :-----------------: | :-------------------------: | :-----------------------: | :--------------------: | :--------------------: |
| [`Contingency`](@ref)     | Contingency table   |         ✓          |             ✓              |            ✓             |           ✓           |           ✓           |
| [`UniqueElements`](@ref)  | Frequencies         |         ✓          |             ✓              |            ✓             |           ✓           |           ✖           |
| [`ValueHistogram`](@ref)  | Binning (histogram) |         ✓          |             ✓              |            ✓             |           ✓           |           ✖           |
| [`OrdinalPatterns`](@ref) | Ordinal patterns    |         ✓          |             ✓              |            ✓             |           ✓           |           ✖           |
| [`Dispersion`](@ref)      | Dispersion patterns |         ✓          |             ✓              |            ✓             |           ✓           |           ✖           |
"""
function mutualinfo(measure::MutualInformation, est::ProbabilitiesEstimator, x, y)
    return estimate(measure, est, x, y)
end

function estimate(est::ProbabilitiesEstimator, x, y)
    estimate(MIShannon(), est, x, y)
end

"""
    mutualinfo([measure::MutualInformation], est::DifferentialEntropyEstimator, x, y)

Estimate the mutual information `measure` between `x` and `y` by a sum of three
entropy terms, without any bias correction, using any [`DifferentialEntropyEstimator`](@ref)
compatible with multivariate data. See examples
[here](@ref example_mi_DifferentialEntropyEstimator). If `measure` is not given, then the
default is `MIShannon()`.

!!! note
    [`DifferentialEntropyEstimator`](@ref)s have their own `base` field which is not
    used here. Instead, this method creates a copy of `est` internally,
    where `est.base` is replaced by `measure.e.base`. Therefore, use `measure` to
    control the "unit" of the mutual information.

## Estimators

Some [`MutualInformation`](@ref) measures can be computed using a [`DifferentialEntropyEstimator`](@ref),
provided it supports multivariate input data. These estimators compute mutual information as a sum of
of entropy terms (with different dimensions), without any bias correction.

| Estimator                        | Principle         | [`MIShannon`](@ref) | [`MITsallisFuruichi`](@ref) | [`MITsallisMartin`](@ref) | [`MIRenyiJizba`](@ref) | [`MIRenyiSurbu`](@ref) |
| -------------------------------- | ----------------- | :-----------------: | :-------------------------: | :-----------------------: | :--------------------: | :--------------------: |
| [`Kraskov`](@ref)                | Nearest neighbors |         ✓          |              x              |             x             |           x            |           x            |
| [`Zhu`](@ref)                    | Nearest neighbors |         ✓          |              x              |             x             |           x            |           x            |
| [`ZhuSingh`](@ref)               | Nearest neighbors |         ✓          |              x              |             x             |           x            |           x            |
| [`Gao`](@ref)                    | Nearest neighbors |         ✓          |              x              |             x             |           x            |           x            |
| [`Goria`](@ref)                  | Nearest neighbors |         ✓          |              x              |             x             |           x            |           x            |
| [`Lord`](@ref)                   | Nearest neighbors |         ✓          |              x              |             x             |           x            |           x            |
| [`LeonenkoProzantoSavani`](@ref) | Nearest neighbors |         ✓          |              x              |             x             |           x            |           x            |
"""
function mutualinfo(est::DifferentialEntropyEstimator, x, y)
    return estimate(est, x, y)
end

# Internal method for compatibility with `independence`
estimate(est::DifferentialEntropyEstimator, x, y) = estimate(MIShannon(), est, x, y)

"""
    mutualinfo([measure::MutualInformation], est::MutualInformationEstimator, x, y)

Estimate the mutual information `measure` between `x` and `y` using the
dedicated [`MutualInformationEstimator`](@ref) `est`.
See examples [here](@ref example_mi_MutualInformationEstimator).
If `measure` is not given, then the default is `MIShannon()`.

## Estimators

Dedicated [`MutualInformationEstimator`](@ref)s are either discrete, continuous,
or a mixture of both. Typically, these estimators apply bias correction.

| Estimator                      |    Type    | [`MIShannon`](@ref) |
| ------------------------------ | :--------: | :-----------------: |
| [`GaussanMI`](@ref)            | Parametric |         ✓          |
| [`KSG1`](@ref)                 | Continuous |         ✓          |
| [`KSG2`](@ref)                 | Continuous |         ✓          |
| [`GaoKannanOhViswanath`](@ref) |   Mixed    |         ✓          |
| [`GaoOhViswanath`](@ref)       | Continuous |         ✓          |
"""
function mutualinfo(measure::MIShannon, est::MutualInformationEstimator, x, y)
    return estimate(MIShannon(measure.e), est, x, y)
end

# Internal method for compatibility with `independence`
estimate(est::MutualInformationEstimator, x, y) = estimate(MIShannon(), est, x, y)

# Generic 3H-formulation of mutual information.
function marginal_entropies_mi3h(measure::MutualInformation, est, x, y)
    e = measure.e
    X = StateSpaceSet(x)
    Y = StateSpaceSet(y)
    XY = StateSpaceSet(X, Y)
    HX = entropy(e, est, X)
    HY = entropy(e, est, Y)
    HXY = entropy(e, est, XY)
    return HX, HY, HXY
end

# Override some definitions, because the estimator behaviour need to be adjusted
# for multiple input variables.
const WellDefinedMIShannonProbEsts{m, D} = Union{
    OrdinalPatterns{m},
    ValueHistogram{<:FixedRectangularBinning{D}},
    ValueHistogram{<:RectangularBinning{T}},
    Dispersion
} where {m, D, T}

function marginal_entropies_mi3h(measure::MutualInformation,
        est::WellDefinedMIShannonProbEsts{m, D}, x, y) where {m, D}
    eX, eY = marginal_encodings(est, x, y)
    eXY = StateSpaceSet(eX, eY)
    e = measure.e
    HX = entropy(e, UniqueElements(), eX)
    HY = entropy(e, UniqueElements(), eY)
    HXY = entropy(e, UniqueElements(), eXY)
    return HX, HY, HXY
end
