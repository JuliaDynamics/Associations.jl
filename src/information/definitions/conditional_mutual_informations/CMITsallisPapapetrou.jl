using ComplexityMeasures: Tsallis
import ComplexityMeasures: log_with_base
export CMITsallisPapapetrou

"""
    CMITsallisPapapetrou <: ConditionalMutualInformation
    CMITsallisPapapetrou(; base = 2, q = 1.5)

The Tsallis-Papapetrou conditional mutual information [Papapetrou2020](@cite).

## Usage

- Use with [`association`](@ref) to compute the raw Tsallis-Papapetrou conditional mutual information
    using of of the estimators listed below.
- Use with [`independence`](@ref) to perform a formal hypothesis test for pairwise conditional 
    independence using the Tsallis-Papapetrou conditional mutual information.

## Compatible estimators

- [`JointProbabilities`](@ref)

## Definition

Tsallis-Papapetrou conditional mutual information is defined as 

```math
I_T^q(X, Y \\mid Z) = \\frac{1}{1 - q} \\left( 1 - \\sum_{XYZ} \\frac{p(x, y, z)^q}{p(x \\mid z)^{q-1} p(y \\mid z)^{q-1} p(z)^{q-1}} \\right).
```
"""
Base.@kwdef struct CMITsallisPapapetrou{B, Q} <: ConditionalMutualInformation
    base::B = 2
    q::Q = 1.5
end

function association(definition::CMITsallisPapapetrou, pxyz::Probabilities{T, 3}) where T
    (; base, q) = definition

    dx, dy, dz = size(pxyz)
    pxz = marginal(pxyz, dims = [1, 3])
    pyz = marginal(pxyz, dims = [2, 3])
    pz = marginal(pxyz, dims = 3)

    cmi = 0.0
    pq = q-1
    for k in 1:dz
        pzₖ = pz[k]
        for j in 1:dy
            pyⱼzₖ = pyz[j, k]
            for i in 1:dx
                pxᵢzₖ = pxz[i, k]
                pxᵢyⱼzₖ = pxyz[i, j, k]
                if pxᵢzₖ != 0.0 && pyⱼzₖ != 0.0 && pzₖ != 0.0
                    cmi +=pxᵢyⱼzₖ / (pxᵢzₖ^pq * pyⱼzₖ^pq * pzₖ^pq)
                end
            end
        end
    end
    return 1 / (1 - q) * (1 - cmi)
end
