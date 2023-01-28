
export ConditionalMutualInformationEstimator
export ConditionalMutualInformation
export condmutualinfo

"""
    ConditionalMutualInformation <: InformationMeasure
    CMI # alias

The supertype of all conditional mutual informations.
"""
abstract type ConditionalMutualInformation{E} <: InformationMeasure end
const CMI{E} = ConditionalMutualInformation{E}

"""
    ConditionalMutualInformationEstimator <: InformationEstimator
    CMIEstimator # alias

The supertype of all conditional mutual information estimators.

## Subtypes

- [`FPVP`](@ref).
- [`PoczosSchneiderCMI`](@ref).
- [`Rahimzamani`](@ref).
- [`MesnerShalisi`](@ref).
"""
abstract type ConditionalMutualInformationEstimator end
const CMIEstimator = ConditionalMutualInformationEstimator

"""
    condmutualinfo([measure::CMI], est::CMIEstimator, x, y, z) → cmi::Real

Estimate a conditional mutual information (CMI) of some kind (specified by `measure`),
between `x` and `y`, given `z`, using the given dedicated
[`ConditionalMutualInformationEstimator`](@ref), which may be discrete, continuous or
mixed.
"""
condmutualinfo(args...; kwargs...) = estimate(args...; kwargs...)

function condmutualinfo(est::ConditionalMutualInformationEstimator, x, y, z)
    return estimate(CMIShannon(), est, x, y, z)
end

include("CMIShannon.jl")
include("CMIRenyiSarbu.jl")
include("CMIRenyiJizba.jl")
include("CMIRenyiPoczos.jl")
include("estimators/estimators.jl")

# Default to Shannon mutual information.
"""
    condmutualinfo([measure::CMI], est::ProbabilitiesEstimator, x, y, z) → cmi::Real ∈ [0, a)

Estimate the conditional mutual information (CMI) `measure` between `x` and `y` given `z`
using a sum of entropy terms, without any bias correction, using the provided
[`ProbabilitiesEstimator`](@ref) `est`.
If `measure` is not given, then the default is `CMIShannon()`.

With a [`ProbabilitiesEstimator`](@ref), the returned `cmi` is guaranteed to be
non-negative.
"""
function condmutualinfo(est::ProbabilitiesEstimator, x, y, z)
    return estimate(CMIShannon(), est, x, y, z)
end

"""
    condmutualinfo([measure::CMI], est::DifferentialEntropyEstimator, x, y, z) → cmi

Estimate the conditional mutual information (CMI) `measure` between `x` and `y` using
a sum of entropy terms, without any bias correction, using the provided
[`DifferentialEntropyEstimator`](@ref) `est` (which must support multivariate data).
If `measure` is not given, then the default is `CMIShannon()`.
"""
function condmutualinfo(est::DifferentialEntropyEstimator, x, y, z)
    return estimate(CMIShannon(), est, x, y, z)
end

"""
    condmutualinfo([measure::CMI], est::MutualInformationEstimator, x, y, z) → cmi::Real

Estimate the conditional mutual information (CMI) `measure` between `x` and `y` using
a difference of mutual information terms, without any bias correction, using the provided
[`MutualInformationEstimator`](@ref) `est`, which may be continuous/differential,
discrete or mixed.
If `measure` is not given, then the default is `CMIShannon()`.
"""
function condmutualinfo(est::MutualInformationEstimator, x, y, z)
    return estimate(CMIShannon(), est, x, y, z)
end

# Generic H4-formulation of CMI
function marginal_entropies_cmi4h(measure::ConditionalMutualInformation, est, x, y, z)
    e = measure.e
    Z = Dataset(z)
    Y = Dataset(y)
    X = Dataset(x)
    XZ = Dataset(X, Z)
    YZ = Dataset(Y, Z)
    XYZ = Dataset(X, Y, Z)

    HXZ = entropy(e, est, XZ)
    HYZ = entropy(e, est,YZ)
    HXYZ = entropy(e, est, XYZ)
    HZ = entropy(e, est, Z)
    return HXZ, HYZ, HXYZ, HZ
end

# Override some definitions, because the estimator behaviour need to be adjusted
# for multiple input variables.
const WellDefinedCMIShannonProbEsts{m, D} = Union{
    SymbolicPermutation{m},
    ValueHistogram{<:FixedRectangularBinning{D}},
    ValueHistogram{<:RectangularBinning{T}},
    Dispersion
} where {m, D, T}

function marginal_entropies_cmi4h(measure::Union{CMIShannon, CMIRenyiSarbu},
        est::WellDefinedCMIShannonProbEsts{m, D},
        x, y, z) where {m, D}
    e = measure.e
    eX, eY, eZ = marginal_encodings(est, x, y, z)
    eXZ = Dataset(eX, eZ)
    eYZ = Dataset(eY, eZ)
    eXYZ = Dataset(eX, eY, eZ)

    HXZ = entropy(e, CountOccurrences(), eXZ)
    HYZ = entropy(e, CountOccurrences(), eYZ)
    HXYZ = entropy(e, CountOccurrences(), eXYZ)
    HZ = entropy(e, CountOccurrences(), eZ)
    return HXZ, HYZ, HXYZ, HZ
end

mi_measure(m::CMIShannon) = MIShannon(m.e)
function estimate(measure::CMI, est::MutualInformationEstimator, x, y, z)
    X = Dataset(x)
    Y = Dataset(y)
    Z = Dataset(z)
    YZ = Dataset(Y, Z)
    m = mi_measure(measure)
    return mutualinfo(m, est, X, YZ) - mutualinfo(m, est, X, Y)
end
