export ConditionalMutualInformation

"""
    CondiitionalMutualInformation

Abstract type for all mutual information measures.

## Concrete implementations

- [`CMIShannon`](@ref)
- [`CMITsallis`](@ref)
- [`CMIRenyiJizba`](@ref)
- [`CMIRenyiSarbu`](@ref)
- [`CMIRenyiPoczos`](@ref)

See also: [`ConditionalMutualInformationEstimator`](@ref)
"""
abstract type ConditionalMutualInformation <: MultivariateInformationMeasure end

min_inputs_vars(::ConditionalMutualInformation) = 3
max_inputs_vars(::ConditionalMutualInformation) = 3


# Generic H4-formulation of CMI
function marginal_entropies_cmi4h_differential(est::EntropyDecomposition{<:ConditionalMutualInformation, <:DifferentialInfoEstimator}, x, y, z)
    Z = StateSpaceSet(z)
    Y = StateSpaceSet(y)
    X = StateSpaceSet(x)
    XZ = StateSpaceSet(X, Z)
    YZ = StateSpaceSet(Y, Z)
    XYZ = StateSpaceSet(X, Y, Z)

    modified_est = estimator_with_overridden_parameters(est.definition, est.est)
    HXZ = entropy(modified_est, XZ)
    HYZ = entropy(modified_est, YZ)
    HXYZ = entropy(modified_est, XYZ)
    HZ = entropy(modified_est, Z)

    return HXZ, HYZ, HXYZ, HZ
end

function marginal_entropies_cmi4h_discrete(est::EntropyDecomposition{<:ConditionalMutualInformation, <:DiscreteInfoEstimator}, x, y, z)
    # Encode marginals to integers based on the outcome space.
    eX, eY, eZ = codified_marginals(est.discretization, x, y, z)
    eXZ = StateSpaceSet(eX, eZ)
    eYZ = StateSpaceSet(eY, eZ)
    eXYZ = StateSpaceSet(eX, eY, eZ)
 
    # The outcome space is no longer relevant from this point on. We're done discretizing, 
    # so now we can just count (i.e. use `UniqueElements` as the outcome space).
    o = UniqueElements()

    modified_est = estimator_with_overridden_parameters(est.definition, est.mest)
    HXZ = information(modified_est, est.pest, o, eXZ)
    HYZ = information(modified_est, est.pest, o, eYZ)
    HXYZ = information(modified_est, est.pest, o, eXYZ)
    HZ = information(modified_est, est.pest, o, eZ)

    return HXZ, HYZ, HXYZ, HZ
end

function information(est::CMIDecomposition{<:ConditionalMutualInformation}, x, y, z)
    return information(est.est, x, y, z)
end

include("CMIShannon.jl")
include("CMITsallis.jl")
include("CMIRenyiJizba.jl")
include("CMIRenyiPoczos.jl")
include("CMIRenyiSarbu.jl")
