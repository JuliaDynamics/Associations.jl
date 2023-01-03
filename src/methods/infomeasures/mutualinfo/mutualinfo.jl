export MutualInformation
export MutualInformationEstimator
export MutualInformationDefinition
export mutualinfo

""" The supertype of all mutual information measures """
abstract type MutualInformation{E, D} <: InformationMeasure end

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

Estimate the mutual information `measure` (either [`MIShannon`](@ref) or
[`MITsallis`](@ref)) between `x` and `y` using the provided estimator `est`.

The online documentation provides an [overview table](@ref mutualinfo_overview) of
compatible measures/definition and estimators.
"""
mutualinfo(args...; kwargs...) = estimate(args...; kwargs...)

include("definitions/definitions.jl")
include("MIShannon.jl")
include("MITsallis.jl")
include("MIRenyi.jl")

include("common_dispatch.jl")

include("estimators/estimators.jl")

# Default to Shannon mutual information.
mutualinfo(est::ProbOrDiffEst, x, y) = estimate(MIShannon(), est, x, y)
