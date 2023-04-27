using LinearAlgebra: eigen, inv
using StateSpaceSets: AbstractStateSpaceSet, StateSpaceSet
using StaticArrays: SUnitRange, SOneTo
using Statistics: mean, std, cor

export cca_cor
function cca_cor(x::AbstractStateSpaceSet{D1}, x̂::AbstractStateSpaceSet{D2}) where {D1, D2}
    # Sample covariance matrices
    data = StateSpaceSet(x, x̂)
    μ, R = fastmean_and_cov(data) # This comes in the Gao2017.jl file and is faster than converting to matrix first.
    Rxx̂ = R[1:D1, (D1+1):(D1+D2)]
    Rx̂x = R[(D1+1):(D1+D2), 1:D1]
    Rxx = R[1:D1, 1:D1]
    Rx̂x̂ = R[(D1+1):(D1+D2), (D1+1):(D1+D2)]
    Rx = inv(Rxx) * Rxx̂ * inv(Rx̂x̂) * Rx̂x
    Rx̂ = inv(Rx̂x̂) * Rx̂x * inv(Rxx) * Rxx̂

    # Eigenvalues (plural) λx and eigenvectors (plural) vsx, and the same for x̂
    λx, vsx = eigen(Rx)
    λx̂, vsx̂ = eigen(Rx̂)

    # Scaling vectors(i.e. weights for linear combinations) are the eigenvectors
    # associated with the largest eigenvalue
    ix = last(findmax(λx))
    ix̂ = last(findmax(λx̂))
    vx = vsx[SOneTo(D1), ix]
    vx̂ = vsx̂[SOneTo(D2), ix̂]

    # Calculates scores `sx` and `sy` (ie. transform data into reduced 2D space)
    # by a linear combination
    sx = [vx ⋅ xᵢ for xᵢ in x] # dot product
    sx̂ = [vx̂ ⋅ x̂ᵢ for x̂ᵢ in x̂] # dot product
    return sx, sx̂
    # We don't need any information beyond the correlation between the reduced variables.
    #  in SX is now a linear combination of the variables in x
    # and each variable in SX̂ is a linear combination of the variables in x̂
    #return cor(sx, sx̂)
    # # Standardize and compute z-scores
    # μsx, σsx = mean(sx), std(sx)
    # μsx̂, σsx̂ = mean(sx̂), std(sx̂)

    # # Standardized scores
    # ssx = (sx .- μsx) / σsx
    # ssx̂ = (sx̂ .- μsx̂) / σsx̂
end
