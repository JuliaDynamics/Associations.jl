
export ConditionalMutualInformationEstimator
export ConditionalMutualInformation
export condmutualinfo

"""
    ConditionalMutualInformation <: AssociationMeasure
    CMI # alias

The supertype of all conditional mutual information measures. Concrete subtypes are

- [`CMIShannon`](@ref)
- [`CMIRenyiJizba`](@ref)
- [`CMIRenyiPoczos`](@ref)
"""
abstract type ConditionalMutualInformation{E} <: MultivariateInformationMeasure end
const CMI{E} = ConditionalMutualInformation{E}

"""
    ConditionalMutualInformationEstimator <: MultivariateInformationMeasureEstimator
    CMIEstimator # alias

The supertype of all conditional mutual information estimators.

## Subtypes

- [`FPVP`](@ref).
- [`PoczosSchneiderCMI`](@ref).
- [`Rahimzamani`](@ref).
- [`MesnerShalizi`](@ref).
"""
abstract type ConditionalMutualInformationEstimator <: MultivariateInformationMeasureEstimator end
const CMIEstimator = ConditionalMutualInformationEstimator

condmutualinfo(args...; kwargs...) = estimate(args...; kwargs...)

const CMI_ESTIMATORS = Union{
    OutcomeSpace,
    DifferentialInfoEstimator,
    MutualInformationEstimator,
    ConditionalMutualInformationEstimator
}
function estimate(measure::CMI, est::CMI_ESTIMATORS, x)
    txt = "`condmutualinfo` takes three input vectors/StateSpaceSets. Only one was given."
    throw(ArgumentError(txt))
end
function estimate(measure::CMI, est::CMI_ESTIMATORS, x, y)
    txt = "`condmutualinfo` takes three input vectors/StateSpaceSets. Only two were given."
    throw(ArgumentError(txt))
end

"""
    condmutualinfo([measure::CMI], est::CMIEstimator, x, y, z) → cmi::Real

Estimate a conditional mutual information (CMI) of some kind (specified by `measure`),
between `x` and `y`, given `z`, using the given dedicated
[`ConditionalMutualInformationEstimator`](@ref), which may be discrete, continuous or
mixed.

## Estimators

| Estimator                    | Principle         | [`CMIShannon`](@ref) | [`CMIRenyiPoczos`](@ref) |
| ---------------------------- | ----------------- | :------------------: | :----------------------: |
| [`FPVP`](@ref)               | Nearest neighbors |          ✓          |            x             |
| [`MesnerShalizi`](@ref)      | Nearest neighbors |          ✓          |            x             |
| [`Rahimzamani`](@ref)        | Nearest neighbors |          ✓          |            x             |
| [`PoczosSchneiderCMI`](@ref) | Nearest neighbors |          x           |            ✓            |
"""
function condmutualinfo(measure::CMI, est::ConditionalMutualInformationEstimator, x, y, z)
    return estimate(measure, est, x, y, z)
end

function estimate(est::ConditionalMutualInformationEstimator, x, y, z)
    return estimate(CMIShannon(), est, x, y, z)
end


include("CMIShannon.jl")
include("CMIRenyiSarbu.jl")
include("CMIRenyiJizba.jl")
include("CMIRenyiPoczos.jl")
include("estimators/estimators.jl")

# Default to Shannon mutual information.
"""
    condmutualinfo([measure::CMI], est::OutcomeSpace, x, y, z) → cmi::Real ∈ [0, a)

Estimate the conditional mutual information (CMI) `measure` between `x` and `y` given `z`
using a sum of entropy terms, without any bias correction, using the provided
[`OutcomeSpace`](@ref) `est`.
If `measure` is not given, then the default is `CMIShannon()`.

With a [`OutcomeSpace`](@ref), the returned `cmi` is guaranteed to be
non-negative.

## Estimators

| Estimator                    | Principle           | [`CMIShannon`](@ref) | [`CMIRenyiSarbu`](@ref) |
| ---------------------------- | ------------------- | :------------------: | :---------------------: |
| [`CountOccurrences`](@ref)   | Frequencies         |          ✓          |           ✓            |
| [`ValueHistogram`](@ref)     | Binning (histogram) |          ✓          |           ✓            |
| [`OrdinalPatterns`](@ref) | Ordinal patterns    |          ✓          |           ✓            |
| [`Dispersion`](@ref)         | Dispersion patterns |          ✓          |           ✓            |
"""
function condmutualinfo(measure::CMI, est::OutcomeSpace, x, y, z)
    return estimate(measure, est, x, y, z)
end

function estimate(est::OutcomeSpace, x, y, z)
    return estimate(CMIShannon(), est, x, y, z)
end

"""
    condmutualinfo([measure::CMI], est::DifferentialInfoEstimator, x, y, z) → cmi::Real

Estimate the mutual information between `x` and `y` conditioned on `z`, using
the differential version of the given conditional mutual information (CMI) `measure`.
The [`DifferentialInfoEstimator`](@ref) `est` must must support multivariate data.
No bias correction is performed. If `measure` is not given, then the default is
`CMIShannon()`.

!!! note
    [`DifferentialInfoEstimator`](@ref)s have their own `base` field which is not
    used here. Instead, this method creates a copy of `est` internally,
    where `est.base` is replaced by `measure.e.base`. Therefore, use `measure` to
    control the "unit" of the mutual information.

## Estimators

| Estimator                        | Principle         | [`CMIShannon`](@ref) |
| -------------------------------- | ----------------- | :------------------: |
| [`Kraskov`](@ref)                | Nearest neighbors |          ✓          |
| [`Zhu`](@ref)                    | Nearest neighbors |          ✓          |
| [`Gao`](@ref)                    | Nearest neighbors |          ✓          |
| [`Goria`](@ref)                  | Nearest neighbors |          ✓          |
| [`Lord`](@ref)                   | Nearest neighbors |          ✓          |
"""
function condmutualinfo(measure::CMI, est::DifferentialInfoEstimator, x, y, z)
    return estimate(measure, est, x, y, z)
end

function estimate(est::DifferentialInfoEstimator, x, y, z)
    return estimate(CMIShannon(), est, x, y, z)
end

"""
    condmutualinfo([measure::CMI], est::MutualInformationEstimator, x, y, z) → cmi::Real

Estimate the mutual information between `x` and `y` conditioned on `z`, using the
given conditional mutual information (CMI) `measure`, computed as a
a difference of mutual information terms (just the chain rule of mutual information)

```math
\\hat{I}(X; Y | Z) = \\hat{I}(X; Y, Z) - \\hat{I}(X; Z).
```

The [`MutualInformationEstimator`](@ref) `est` may be continuous/differential,
discrete or mixed. No bias correction in performed, except the bias correction
that occurs for each individual mutual information term.
If `measure` is not given, then the default is `CMIShannon()`.

## Estimators

| Estimator                              |    Type    |     Principle     | [`CMIShannon`](@ref) |
| -------------------------------------- | :--------: | :---------------: | :------------------: |
| [`KraskovStögbauerGrassberger1`](@ref) | Continuous | Nearest neighbors |          ✓          |
| [`KraskovStögbauerGrassberger2`](@ref) | Continuous | Nearest neighbors |          ✓          |
| [`GaoKannanOhViswanath`](@ref)         |   Mixed    | Nearest neighbors |          ✓          |
| [`GaoOhViswanath`](@ref)               | Continuous | Nearest neighbors |          ✓          |
"""
function condmutualinfo(measure::CMI, est::MutualInformationEstimator, x, y, z)
    return estimate(measure, est, x, y, z)
end

mi_measure(m::CMIShannon) = MIShannon(m.e)
# Internal methods for `independence`
function estimate(measure::CMI, est::MutualInformationEstimator, x, y, z)
    X = StateSpaceSet(x)
    Y = StateSpaceSet(y)
    Z = StateSpaceSet(z)
    YZ = StateSpaceSet(Y, Z)
    m = mi_measure(measure)
    return mutualinfo(m, est, X, YZ) - mutualinfo(m, est, X, Z)
end

function estimate(est::MutualInformationEstimator, x, y, z)
    return estimate(CMIShannon(), est, x, y, z)
end


# Generic H4-formulation of CMI
function marginal_entropies_cmi4h(measure::ConditionalMutualInformation, est, x, y, z)
    e = measure.e
    Z = StateSpaceSet(z)
    Y = StateSpaceSet(y)
    X = StateSpaceSet(x)
    XZ = StateSpaceSet(X, Z)
    YZ = StateSpaceSet(Y, Z)
    XYZ = StateSpaceSet(X, Y, Z)

    HXZ = information(e, est, XZ)
    HYZ = information(e, est,YZ)
    HXYZ = information(e, est, XYZ)
    HZ = information(e, est, Z)
    return HXZ, HYZ, HXYZ, HZ
end

# Override some definitions, because the estimator behaviour need to be adjusted
# for multiple input variables.
const WellDefinedCMIShannonProbEsts{m, D} = Union{
    OrdinalPatterns{m},
    ValueHistogram{<:FixedRectangularBinning{D}},
    ValueHistogram{<:RectangularBinning{T}},
    Dispersion
} where {m, D, T}

function marginal_entropies_cmi4h(measure::Union{CMIShannon, CMIRenyiSarbu},
        est::WellDefinedCMIShannonProbEsts{m, D},
        x, y, z) where {m, D}
    e = measure.e
    eX, eY, eZ = marginal_encodings(est, x, y, z)
    eXZ = StateSpaceSet(eX, eZ)
    eYZ = StateSpaceSet(eY, eZ)
    eXYZ = StateSpaceSet(eX, eY, eZ)

    HXZ = information(e, CountOccurrences(), eXZ)
    HYZ = information(e, CountOccurrences(), eYZ)
    HXYZ = information(e, CountOccurrences(), eXYZ)
    HZ = information(e, CountOccurrences(), eZ)
    return HXZ, HYZ, HXYZ, HZ
end
