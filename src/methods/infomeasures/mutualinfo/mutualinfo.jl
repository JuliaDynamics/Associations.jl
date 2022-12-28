export MutualInformation
export MutualInformationEstimator
export MutualInformationDefinition
export mutualinfo

""" The supertype of all mutual information like measures """
abstract type MutualInformation <: InformationMeasure end

""" The supertype for mutual information definitions. """
abstract type MutualInformationDefinition <: Definition end

""" The supertype of all dedicated mutual information estimators """
abstract type MutualInformationEstimator <: InformationMeasureEstimator end

"""
    mutualinfo(measure, est::MutualInformationEstimator, x, y)
    mutualinfo([definition], measure, est::DifferentialEntropyEstimator, x, y)
    mutualinfo([definition], measure, est::ProbabilitiesEstimator, x, y)

Estimate  `measure`, a generalized mutual information (MI),
between `x` and `y` using the provided estimator.

- If `est` is a [`MutualInformationEstimator`](@ref), then mutual information is computed
    using some specialized procedure.
- If `est` is an [`DifferentialEntropyEstimator`](@ref) or a [`ProbabilitiesEstimator`], then mutual
    information is computed using the formula specified by `definition`, which is a
    [`MutualInformationDefinition`](@ref) (listed below).

## Supported measures (definitions)

- [`MIShannon`](@ref) (defaults to [`ShannonH3`])
- [`MITsallis`](@ref) (defaults to [`TsallisH3`]).

## Supported estimators (continuous)

Dedicated estimators:

- [`KraskovStögbauerGrassberger1`](@ref), or [`KSG1`](@ref)
- [`KraskovStögbauerGrassberger2`](@ref), or [`KSG2`](@ref)
- [`GaoKannanOhViswanath`](@ref).
- [`GaoOhViswanath`](@ref).

We also support the following [`DifferentialEntropyEstimator`](@ref)s

- [`Kraskov`](@ref)
- [`KozachenkoLeonenko`](@ref)
- [`GaoNaive`](@ref) and [`GaoNaiveCorrected`](@ref)
- [`Goria`](@ref)
- [`Zhu`](@ref)
- [`ZhuSingh`](@ref)
- [`LeonenkoProzantoSavani`](@ref)
- [`Lord`](@ref)

## Supported estimators (discrete)

In principle, any [`ProbabilitiesEstimator`](@ref) that operates on multivariate data
may be used, but you shuld make sure to use one that uses the same discretization for
both the joint and marginal spaces. This is automatically handled if you use:

- **[`ValueHistogram`](@ref)**. Bin visitation frequencies are counted in the joint space
    `XY`, then marginal probabilities are obtained from the joint bin visits.
"""
mutualinfo(args...; kwargs...) = estimate(args...; kwargs...)

include("shannon/MIShannon.jl")
include("tsallis/MITsallis.jl")
include("estimators/estimators.jl")


# estimate(def::MutualInformationDefinition, measure::MutualInformation, est, x, y) =
#     estimate(def, measure, est, x, y)
# estimate(measure::MutualInformation, est, x, y) =
#     estimate(measure, est, x, y)

# # Default to Shannon mutual information
# estimate(est::ProbabilitiesEstimator, x, y) =
#     estimate(MIShannon(), est, x, y)
