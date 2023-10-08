abstract type MutualInformation <: BivariateInformationMeasure end

# Generic 3H-formulation of Shannon mutual information (a scaled sum of these entropies are also used by 
# some of the other mutual information measures, so we define it generically here).
function marginal_entropies_mi3h_differential(est::DifferentialDecomposition{<:MutualInformation, <:DifferentialInfoEstimator}, x, y)
    X = StateSpaceSet(x)
    Y = StateSpaceSet(y)
    XY = StateSpaceSet(X, Y)
    HX = information(est.est, X) # estimates entropy in the X marginal
    HY = information(est.est, Y) # estimates entropy in the Y marginal
    HXY = information(est.est, XY) # estimates entropy in the joint XY space
    return HX, HY, HXY
end

function marginal_entropies_mi3h_discrete(est::DiscreteDecomposition{<:MutualInformation, <:DiscreteInfoEstimator}, x, y)
    # Encode marginals to integers based on the outcome space.
    eX, eY = codified_marginals(est.discretization, x, y)
    eXY = StateSpaceSet(eX, eY)

    # The outcome space is no longer relevant from this point on. We're done discretizing, 
    # so now we can just count (i.e. use `UniqueElements`).
    o = UniqueElements()
    HX = information(est.mest, est.pest, o, eX) # estimates entropy in the X marginal
    HY = information(est.mest, est.pest, o, eY) # estimates entropy in the Y marginal
    HXY = information(est.mest, est.pest, o, eXY) # estimates entropy in the joint XY space
    
    return HX, HY, HXY
end


include("MIShannon.jl")
include("MITsallisMartin.jl")
include("MITsallisFuruichi.jl")
include("MIRenyiSarbu.jl")
include("MIRenyiJizba.jl")