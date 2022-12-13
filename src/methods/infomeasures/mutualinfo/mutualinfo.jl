export MutualInformation
export MutualInformationEstimator
export mutualinfo

""" The supertype of all mutual information like measures """
abstract type MutualInformation <: InformationMeasure end

""" The supertype of all dedicated mutual information estimators """
abstract type MutualInformationEstimator <: InformationMeasureEstimator end


"""
    mutualinfo(definition::MIShannon, est::ProbabilitiesEstimator, x, y) → MI::Real
    mutualinfo(definition::MITsallis, est::ProbabilitiesEstimator, x, y) → MI::Real

    mutualinfo(definition::MIShannonDifferential, est::EntropyEstimator, x, y) → MI::Real
    mutualinfo(definition::MIShannonDifferential, est::MutualInformationEstimator, x, y) → MI::Real

Estimate the generalized mutual information (MI) (using the formula and logarithm base
specified by `definition`) between `x` and `y`, using the provided estimator.

The first set of signatures is for discrete MI. The second set of signatures is for
differential MI. For a full list of compatible definitions and estimators, see the
online documentation.

Returns `MI`, the mutual information estimate, whose interpretation depends on the
combination of `definition` and `est`.

## Supported definitions

Generalized mutual information-like quantities are abundant in the literature.
Sometimes, different authors give different definitions of divergence measures with the
same name. We currently support the following measures (some of which may be tweaked
according to multiple definitions):

- [`MIShannon`](@ref). Discrete Shannon mutual information.
- [`MITsallis`](@ref). Discrete Tsallis mutual information.
- [`MIShannonDifferential`](@ref). Differential Shannon mutual information.
"""
mutualinfo(def, est, x, y) = estimate(def, est, x, y)

include("shannon/MIShannon.jl")
include("shannon/MIShannonDifferential.jl")
include("tsallis/MITsallis.jl")

include("estimators/estimators.jl")
