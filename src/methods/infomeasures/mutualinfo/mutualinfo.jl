export MutualInformation
export MutualInformationEstimator
export MutualInformationDefinition
export mutualinfo

""" The supertype of all mutual information measures """
abstract type MutualInformation{E} <: InformationMeasure end

""" The supertype of all dedicated mutual information estimators """
abstract type MutualInformationEstimator <: InformationMeasureEstimator end

# There are many ways of defining mutual information. Moreover, the definitions
# differ for different types of base `EntropyDefinition`s. Therefore, we dispatch
# on subtypes of `MutualInformationDefinition`.
""" The supertype for mutual information definitions. """
abstract type MutualInformationDefinition <: Definition end

""" The supertype for all H3-type (three entropies) decomposition of mutual information. """
abstract type MIH3 <: MutualInformationDefinition end

"""
    mutualinfo(measure::MutualInformation, est::MutualInformationEstimator, x, y)
    mutualinfo(measure::MutualInformation, est::DifferentialEntropyEstimator, x, y)
    mutualinfo(measure::MutualInformation, est::ProbabilitiesEstimator, x, y)
    mutualinfo(measure::MutualInformation, c::ContingencyMatrix)

Estimate the mutual information `measure` (either [`MIShannon`](@ref) or
[`MITsallis`](@ref), ) between `x` and `y` using the provided estimator `est`.
Alternatively, compute mutual information from a pre-computed [`ContingencyMatrix`](@ref).

Compatible measures/definitions and estimators are listed in the
[online documentation](@ref mutualinfo_overview).
"""
mutualinfo(args...; kwargs...) = estimate(args...; kwargs...)

include("MIShannon.jl")
include("MITsallisFuruichi.jl")
include("MITsallisMartin.jl")
include("MIRenyiSarbu.jl")
include("MIRenyiJizba.jl")

include("estimators/estimators.jl")

# Default to Shannon mutual information.

"""
    mutualinfo(measure::MutualInformation], est::ContingencyMatrix)

Estimate the discrete version of the given [`MutualInformation`](@ref) `measure` from
its direct definition (double-sum), using the probabilities from a pre-computed
[`ContingencyMatrix`](@ref).

See the [online documentation](@ref contingency_matrix_mi) for a list of
compatible measures.
"""
mutualinfo(c::ContingencyMatrix) = estimate(MIShannon(), c)

"""
    mutualinfo([measure::MutualInformation], est::ProbabilitiesEstimator, x, y)

Estimate the mutual information `measure` between `x` and `y` by a sum of three
entropy terms, without any bias correction, using the provided
[`ProbabilitiesEstimator`](@ref) `est`. If `measure` is not given, then the default
is `MIShannon()`.

Joint and marginal probabilities are computed by jointly discretizing `x` and `y` using
the approach given by `est`, and obtaining marginal distributions from the joint
distribution.

This only works for estimators that have an implementation for
[`marginal_encodings`](@ref). See the
[online documentation](@ref dedicated_probabilities_estimators_mi) for a list of
compatible measures.
"""
mutualinfo(est::ProbabilitiesEstimator, x, y) = estimate(MIShannon(), est, x, y)

"""
    mutualinfo([measure::MutualInformation], est::DifferentialEntropyEstimator, x, y)

Estimate the mutual information `measure` between `x` and `y` by a sum of three
entropy terms, without any bias correction, using any [`DifferentialEntropyEstimator`](@ref)
compatible with multivariate data. If `measure` is not given, then the default
is `MIShannon()`.

See the [online documentation](@ref dedicated_diffentropy_estimators_mi) for a list of
compatible measures.
"""
mutualinfo(est::DifferentialEntropyEstimator, x, y) = estimate(MIShannon(), est, x, y)

"""
    mutualinfo(measure::MutualInformation, est::MutualInformationEstimator, x, y)

Estimate the mutual information `measure` between `x` and `y` using the dedicated
[`MutualInformationEstimator`](@ref) `est`, which can be either discrete, continuous,
or a mixture of both, and typically involve some bias correction.
If `measure` is not given, then the default is `MIShannon()`.

See the [online documentation](@ref dedicated_estimators_mi) for a list of
compatible measures.
"""
mutualinfo(est::MutualInformationEstimator, x, y) = estimate(MIShannon(), est, x, y)

# Generic 3H-formulation of mutual information.
function marginal_entropies_mi3h(measure::MutualInformation, est, x, y)
    e = measure.e
    X = Dataset(x)
    Y = Dataset(y)
    XY = Dataset(X, Y)
    HX = entropy(e, est, X)
    HY = entropy(e, est, Y)
    HXY = entropy(e, est, XY)
    return HX, HY, HXY
end

# Override some definitions, because the estimator behaviour need to be adjusted
# for multiple input variables.
const WellDefinedMIShannonProbEsts{m, D} = Union{
    SymbolicPermutation{m},
    ValueHistogram{<:FixedRectangularBinning{D}},
    ValueHistogram{<:RectangularBinning{T}},
    Dispersion
} where {m, D, T}

function marginal_entropies_mi3h(measure::MutualInformation,
        est::WellDefinedMIShannonProbEsts{m, D}, x, y) where {m, D}
    eX, eY = marginal_encodings(est, x, y)
    eXY = Dataset(eX, eY)
    e = measure.e
    HX = entropy(e, CountOccurrences(), eX)
    HY = entropy(e, CountOccurrences(), eY)
    HXY = entropy(e, CountOccurrences(), eXY)
    return HX, HY, HXY
end
