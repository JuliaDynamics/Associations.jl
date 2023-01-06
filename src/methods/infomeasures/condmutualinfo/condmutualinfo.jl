
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

- [`FrenzelPompeVelmejkaPalus`](@ref).
- [`PoczosSchneiderCMI`](@ref).
- [`Rahimzamani`](@ref).
- [`MesnerShalisi`](@ref).
"""
abstract type ConditionalMutualInformationEstimator end
const CMIEstimator = ConditionalMutualInformationEstimator

"""
    condmutualinfo(measure::CMI, est::CMIEstimator, x, y, z) → cmi::Real
    condmutualinfo(measure::CMI, est::DifferentialEntropyEstimator, x, y, z) → cmi::Real
    condmutualinfo(measure::CMI, est::ProbabilitiesEstimator, x, y, z) → cmi::Real

Estimate a conditional mutual information (CMI) of some kind (specified by `measure`),
between `x` and `y`, given `z`, using the given estimator.

## Definition

CMIs appear in many forms in the scientific literature. We support the following CMIs:

- **[`CMIShannon`](@ref)**: Shannon CMI.
- **[`CMIRenyi`](@ref)**: Renyi CMI.
"""
condmutualinfo(args...; kwargs...) = estimate(args...; kwargs...)

include("CMIShannon.jl")
include("CMIRenyi.jl")
include("estimators/estimators.jl")

# Default to Shannon mutual information.
condmutualinfo(est::ProbOrDiffEst, x, y, z) = estimate(CMIShannon(), est, x, y, z)
condmutualinfo(est::MutualInformationEstimator, x, y, z) =
    estimate(CMIShannon(), est, x, y, z)

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
    Dispersion
} where {m, D}

function marginal_entropies_cmi4h(measure::CMIShannon,
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
