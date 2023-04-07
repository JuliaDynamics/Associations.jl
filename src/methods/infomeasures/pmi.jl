"""
    PMI <: AssociationMeasure

The partial mutual information (PMI) measure of association.
"""
Base.@kwdef struct PMI <: AssociationMeasure
    base::Real = 2
end

export PMI

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
                @show i, j, k, pmi
                if pzₖ > 0.0
                    part1 = (pxyzᵢⱼₖ / pzₖ) / ((pxzᵢₖ / pzₖ ) * (pyzⱼₖ / pzₖ))
                    part2 = (pxzᵢₖ / pzₖ) * sum(pyzⱼₖ > 0 ? py[j] * (pxyzᵢⱼₖ / pyzⱼₖ) : 0 for j = 1:dy)
                    part3 = (pyzⱼₖ / pzₖ) * sum(pxzᵢₖ > 0 ? px[i] * (pxyzᵢⱼₖ / pxzᵢₖ) : 0 for i = 1:dx)
                else
                    part1 = 0.0
                    part2 = 0.0
                    part3 = 0.0
                end
                if part1 > 0 && part2 > 0 && part3 > 0
                    pmi += pxyzᵢⱼₖ * (logb(part1) + logb(part2) + logb(part3))
                end
            end
        end
    end
    return pmi
end
