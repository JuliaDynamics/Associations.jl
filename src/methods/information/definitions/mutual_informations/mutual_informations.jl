export MutualInformation

"""
    MutualInformation

Abstract type for all mutual information measures.

## Concrete implementations

- [`MIShannon`](@ref)
- [`MITsallisMartin`](@ref)
- [`MITsallisFuruichi`](@ref)
- [`MIRenyiJizba`](@ref)
- [`MIRenyiSarbu`](@ref)

See also: [`MutualInformationEstimator`](@ref)
"""
abstract type MutualInformation <: BivariateInformationMeasure end

# Generic 3H-formulation of Shannon mutual information (a scaled sum of these entropies are also used by 
# some of the other mutual information measures, so we define it generically here).
function marginal_entropies_mi3h_differential(est::EntropyDecomposition{<:MutualInformation, <:DifferentialInfoEstimator}, x, y)
    X = StateSpaceSet(x)
    Y = StateSpaceSet(y)
    XY = StateSpaceSet(X, Y)
    modified_est = estimator_with_overridden_parameters(est.definition, est.est)
    HX = information(modified_est, X) # estimates entropy in the X marginal
    HY = information(modified_est, Y) # estimates entropy in the Y marginal
    HXY = information(modified_est, XY) # estimates entropy in the joint XY space
    return HX, HY, HXY
end

function marginal_entropies_mi3h_discrete(est::EntropyDecomposition{<:MutualInformation, <:DiscreteInfoEstimator}, x, y)
    # Encode marginals to integers based on the outcome space.
    eX::StateSpaceSet, eY::StateSpaceSet = StateSpaceSet.(codified_marginals(est.discretization, x, y))
    eXY::StateSpaceSet = StateSpaceSet(eX, eY)

    # The outcome space is no longer relevant from this point on. We're done discretizing, 
    # so now we can just count (i.e. use `UniqueElements` as the outcome space).
    o = UniqueElements()
    modified_est = estimator_with_overridden_parameters(est.definition, est.est)
    HX = information(modified_est, est.pest, o, eX) # estimates entropy in the X marginal
    HY = information(modified_est, est.pest, o, eY) # estimates entropy in the Y marginal
    HXY = information(modified_est, est.pest, o, eXY) # estimates entropy in the joint XY space
    
    return HX, HY, HXY
end


include("MIShannon.jl")
include("MITsallisMartin.jl")
include("MITsallisFuruichi.jl")
include("MIRenyiSarbu.jl")
include("MIRenyiJizba.jl")