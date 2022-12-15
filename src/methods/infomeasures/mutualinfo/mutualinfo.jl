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

include("measures/MIShannon.jl")
include("definitions/ShannonH3.jl")

include("measures/MITsallis.jl")
include("definitions/TsallisH3.jl")

include("estimators/estimators.jl")

"""
    mutualinfo(measure, est::MutualInformationEstimator, x, y)
    mutualinfo([definition], measure, est::EntropyEstimator, x, y)
    mutualinfo([definition], measure, est::ProbabilitiesEstimator, x, y)

Estimate  `measure`, a generalized mutual information (MI),
between `x` and `y` using the provided estimator.

- If `est` is a [`MutualInformationEstimator`](@ref), then mutual information is computed
    using some specialized procedure.
- If `est` is an [`EntropyEstimator`](@ref) or a [`ProbabilitiesEstimator`], then mutual
    information is computed using some [`MutualInformationDefinition`](@ref) (listed below).

## Supported measures (definitions)

- [`MIShannon`](@ref) (defaults to [`ShannonH3`] definition)
- [`MITsallis`](@ref) (defaults to [`TsallisH3`] definition).

## Supported estimators (continuous)

Dedicated estimators:

- [`KraskovStögbauerGrassberger1`](@ref), or [`KSG1`](@ref)
- [`KraskovStögbauerGrassberger2`](@ref), or [`KSG2`](@ref)
- [`GaoKannanOhViswanath`](@ref).
- [`GaoOhViswanath`](@ref).

We also support the following [`EntropyEstimator`](@ref)s

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
function mutualinfo end
mutualinfo(def::MutualInformationDefinition, measure::MutualInformation, est, x, y) =
    estimate(def, measure, est, x, y)
mutualinfo(measure::MutualInformation, est, x, y) =
    estimate(measure, est, x, y)
