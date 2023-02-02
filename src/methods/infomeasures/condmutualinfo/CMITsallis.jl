"""
    CMITsallis <: ConditionalMutualInformation
    CMITsallis(; q = 1.5, base = 2)
"""
struct CMITsallis{E <: Tsallis} <: ConditionalMutualInformation{E}
    e::E
    function CMITsallis(; base::T = 2, q = 1.5) where {T <: Real}
        e = Tsallis(; base, q)
        new{typeof(e)}(e)
    end
end

function estimate(measure::CMITsallis, est::Contingency{<:ProbabilitiesEstimator}, x...)
    return estimate(measure, contingency_matrix(est.est, x...))
end

function estimate(measure::CMITsallis, est::Contingency{<:Nothing}, x...)
    return estimate(measure, contingency_matrix(x...))
end

function estimate(
        measure::CMITsallis,
        pxyz::ContingencyMatrix{T, 3}) where T
    e = measure.e
    dx, dy, dz = size(pxyz)
    pxz = probabilities(pxyz, dims = [1, 3])
    pyz = probabilities(pxyz, dims = [2, 3])
    pz = probabilities(pxyz, dims = 3)
    cmi = 0.0
    log0 = log_with_base(e.base)
    for k in 1:dz
        pzₖ = pz[k]
        for j in 1:dy
            pyⱼzₖ = pyz[j, k]
            for i in 1:dx
                pxᵢzₖ = pxz[i, k]
                pxᵢyⱼzₖ = pxyz[i, j, k]
                if pxᵢyⱼzₖ != 0.0
                    cmi += pxᵢyⱼzₖ * log0((pzₖ * pxᵢyⱼzₖ) / (pxᵢzₖ * pyⱼzₖ))
                end
            end
        end
    end
    return cmi
end
