abstract type ConditionalMutualInformation <: MultivariateInformationMeasure end

min_inputs_vars(::ConditionalMutualInformation) = 3
max_inputs_vars(::ConditionalMutualInformation) = 3


# Generic H4-formulation of CMI
function marginal_entropies_cmi4h_differential(est::DifferentialDecomposition{<:ConditionalMutualInformation, <:DifferentialInfoEstimator}, x, y, z)
    Z = StateSpaceSet(z)
    Y = StateSpaceSet(y)
    X = StateSpaceSet(x)
    XZ = StateSpaceSet(X, Z)
    YZ = StateSpaceSet(Y, Z)
    XYZ = StateSpaceSet(X, Y, Z)

    HXZ = entropy(est.est, XZ)
    HYZ = entropy(est.est, YZ)
    HXYZ = entropy(est.est, XYZ)
    HZ = entropy(est.est, Z)

    return HXZ, HYZ, HXYZ, HZ
end

function marginal_entropies_cmi4h_discrete(est::DifferentialDecomposition{<:ConditionalMutualInformation, <:DifferentialInfoEstimator}, x, y, z)
    # Encode marginals to integers based on the outcome space.
    eX, eY, eZ = codified_marginals(est, x, y, z)
    eXZ = StateSpaceSet(eX, eZ)
    eYZ = StateSpaceSet(eY, eZ)
    eXYZ = StateSpaceSet(eX, eY, eZ)
 
    # The outcome space is no longer relevant from this point on. We're done discretizing, 
    # so now we can just count (i.e. use `UniqueElements`).
    o = UniqueElements()
    HXZ = information(est.mest, est.pest, o, eXZ)
    HYZ = information(est.mest, est.pest, o, eYZ)
    HXYZ = information(est.mest, est.pest, o, eXYZ)
    HZ = information(est.mest, est.pest, o, eZ)

    return HXZ, HYZ, HXYZ, HZ
end


include("CMIShannon.jl")
include("CMITsallis.jl")
include("CMIRenyiJizba.jl")
include("CMIRenyiPoczos.jl")
include("CMIRenyiSarbu.jl")