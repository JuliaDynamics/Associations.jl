export CMIShannon
import ComplexityMeasures: log_with_base

"""
    CMIShannon <: ConditionalMutualInformation
    CMIShannon(; base = 2)

The Shannon conditional mutual information (CMI) ``I^S(X; Y | Z)``.

## Supported definitions

Consider random variables ``X \\in \\mathbb{R}^{d_X}`` and
``Y \\in \\mathbb{R}^{d_Y}``, given ``Z \\in \\mathbb{R}^{d_Z}``. The Shannon
conditional mutual information is defined as

```math
\\begin{align*}
I(X; Y | Z)
&= H^S(X, Z) + H^S(Y, z) - H^S(X, Y, Z) - H^S(Z) \\\\
&= I^S(X; Y, Z) + I^S(X; Y)
\\end{align*},
```

where ``I^S(\\cdot; \\cdot)`` is the Shannon mutual information [`MIShannon`](@ref),
and ``H^S(\\cdot)`` is the [`Shannon`](@ref) entropy.

Differential Shannon CMI is obtained by replacing the entropies by
differential entropies.

See also: [`condmutualinfo`](@ref).
"""
struct CMIShannon{E <: Shannon} <: ConditionalMutualInformation{E}
    e::E
    function CMIShannon(; base::T = 2) where {T <: Real}
        e = Shannon(; base)
        new{typeof(e)}(e)
    end
end

function estimate(measure::CMIShannon, est::ProbOrDiffEst, x, y, z)
    e = measure.e
    X = Dataset(x)
    Y = Dataset(y)
    Z = Dataset(z)
    XZ = Dataset(X, Z)
    YZ = Dataset(Y, Z)
    XYZ = Dataset(X, Y, Z)
    cmi = entropy(e, est, XZ) +
        entropy(e, est, YZ)  -
        entropy(e, est, XYZ) -
        entropy(e, est, Z)
    return cmi
end


function estimate(
        measure::CMIShannon,
        pxyz::ContingencyMatrix{T, 3}) where T
    e = measure.e
    dx, dy, dz = size(pxyz)
    pxz = dropdims(sum(pxyz, dims = 2), dims = 2)
    pyz = dropdims(sum(pxyz, dims = 1), dims = 1)
    pz = probabilities(pxyz, 3)
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
