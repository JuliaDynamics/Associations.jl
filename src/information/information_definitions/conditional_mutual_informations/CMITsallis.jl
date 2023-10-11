using ComplexityMeasures: Tsallis
import ComplexityMeasures: log_with_base

export CMITsallis

"""
    CMITsallis <: ConditionalMutualInformation
    CMITsallis(; base = 2, q = 1.5)
"""
Base.@kwdef struct CMITsallis{B, Q} <: ConditionalMutualInformation
    base::B = 2
    q::Q = 1.5 # Todo: check formula. Where does `q` appear?
end

# TODO: this has to be wrong. See the following paper for a definition:
# https://www.sciencedirect.com/science/article/pii/S0378437119317030?casa_token=OjBFqwe6_1UAAAAA:s32zgg706dkO5P8Wuy2x2fWX0WjR09IguDvA8xGt6OTthjGkTAAutiZSiqI49ScyA6nTwomXEw
function information(definition::CMITsallis, pxyz::Probabilities{T, 3}) where T
    (; base, q) = definition

    dx, dy, dz = size(pxyz)
    pxz = marginal(pxyz, dims = [1, 3])
    pyz = marginal(pxyz, dims = [2, 3])
    pz = marginal(pxyz, dims = 3)

    cmi = 0.0
    log0 = log_with_base(base)
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
