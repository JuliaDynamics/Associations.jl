export PMI
export pmi

"""
    PMI <: AssociationMeasure

The partial mutual information (PMI) measure of association.

## Estimation

PMI can be estimated using any estimator that implements [`marginal_encodings`](@ref).

## Properties

For the discrete case, the following identities hold in theory (when estimating PMI, they
may not).

- `PMI(X, Y, Z) >= CMI(X, Y, Z)` (where CMI is the Shannon CMI). Holds in theory, but
    when estimating PMI, the identity may not hold.
- `PMI(X, Y, Z) >= 0`. Holds both in theory and for estimation using
    [`ProbabilitiesEstimator`](@ref)s.
- `X ⫫ Y | Z => PMI(X, Y, Z) = CMI(X, Y, Z) = 0` (in theory, but not necessarily for
    estimation).
"""
Base.@kwdef struct PMI <: AssociationMeasure
    base::Real = 2
end

function pmi(x...)
    return estimate(PMI(), x...)
end

function estimate(measure::PMI, est::Contingency{<:ProbabilitiesEstimator}, x...)
    return estimate(measure, contingency_matrix(est.est, x...))
end

function estimate(measure::PMI, est::Contingency{<:Nothing}, x...)
    return estimate(measure, contingency_matrix(x...))
end

# p and q must have the same outcome space.
function extended_kl_div(p, q)

end

function estimate(
        measure::PMI,
        pxyz::ContingencyMatrix{T, 3}) where T
    # The sums go over *states*, so these are what we iterate over.
    dx, dy, dz = size(pxyz)
    pyxz = probabilities(pxyz, dims = [3, 1, 2])
    px = probabilities(pxyz, dims = [1])
    py = probabilities(pxyz, dims = [2])
    pz = probabilities(pxyz, dims = [2])
    pxz = probabilities(pxyz, dims = [1, 3])
    pxy = probabilities(pxyz, dims = [1, 2])
    pyz = probabilities(pxyz, dims = [2, 3])
    pz = probabilities(pxyz, dims = 3)

    pmi = 0.0
    logb = log_with_base(measure.base)
    for i in 1:dx
        for j in 1:dy
            for k in 1:dz
                pzₖ = pz[k]
                pyzⱼₖ = pyz[j, k]
                pxzᵢₖ = pxz[i, k]
                pxyzᵢⱼₖ = pxyz[i, j, k]
                sy = sum(pyzⱼₖ > 0 ? py[j] * (pxyz[i, j, k]/ pyz[j, k]) : 0 for j = 1:dy)
                sx = sum(pxzᵢₖ > 0 ? px[i] * (pxyz[i, j, k] / pxz[i, k]) : 0 for i = 1:dx)

                sxy = sy * sx
                pxy_z = pxyzᵢⱼₖ / pzₖ
                if sxy > 0 && pxy_z > 0
                    pmi += pxyzᵢⱼₖ * logb(pxy_z / (sy * sx))
                end
            end
        end
    end
    return pmi
end
