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

include("estimators/estimators.jl")

# Default to Shannon mutual information.
mutualinfo(est::ProbOrDiffEst, x, y) = estimate(MIShannon(), est, x, y)

# Generic 3H-formulation of mutual information.
function marginal_entropies_mi3h(measure::MutualInformation, est, x, y)
    e = measure.e
    X = Dataset(x)
    Y = Dataset(y)
    XY = Dataset(x, y)
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
    Dispersion
} where {m, D}

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
