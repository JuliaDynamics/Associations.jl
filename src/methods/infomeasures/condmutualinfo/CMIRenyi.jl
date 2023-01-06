export CMIRenyi

"""
    CMIRenyi <: ConditionalMutualInformation
    CMIRenyi(; base = 2, definition = CMIRenyiSarbu())

Rényi conditional mutual information.

## Discrete dDescription

The discrete definition comes from Sarbu (2014)[^Sarbu2014]).

Assume we observe three discrete random variables ``X``, ``Y`` and ``Z``.
Sarbu (2014) defines discrete conditional Rényi mutual information as the conditional
Rényi ``\\alpha``-divergence between the conditional joint probability mass function
``p(x, y | z)`` and the product of the conditional marginals, ``p(x |z) \\cdot p(y|z)``:

```math
I(X, Y; Z)^R_q =
\\dfrac{1}{q-1} \\sum_{z \\in Z} p(Z = z)
\\log \\left(
    \\sum{x \\in X}\\sum{y \\in Y}
    \\dfrac{p(x, y|z)^q}{\\left( p(x|z)\\cdot p(y|z) \\right)^{q-1}}
\\right)
```

[^Sarbu2014]: Sarbu, S. (2014, May). Rényi information transfer: Partial Rényi transfer
    entropy and partial Rényi mutual information. In 2014 IEEE International Conference
    on Acoustics, Speech and Signal Processing (ICASSP) (pp. 5666-5670). IEEE.

See also: [`condmutualinfo`](@ref).
"""
struct CMIRenyi{E <: Renyi} <: ConditionalMutualInformation{E}
    e::E
    function CMIRenyi(; base = 2, q = 1.5)
        e = Renyi(; base, q)
        new{typeof(e)}(e)
    end
end

function estimate(
        measure::CMIRenyi,
        pxyz::ContingencyMatrix{T, 3}) where T
    e = measure.e
    dx, dy, dz = size(pxyz)
    pxz = dropdims(sum(pxyz, dims = 2), dims = 2)
    pyz = dropdims(sum(pxyz, dims = 1), dims = 1)
    pz = probabilities(pxyz, 3)
    cmi = 0.0
    logb = log_with_base(e.base)
    for k in 1:dz
        pzₖ = pz[k]
        inner = 0.0
        if pzₖ != 0.0
            for j in 1:dy
                pyⱼzₖ = pyz[j, k]
                for i in 1:dx
                    pxᵢzₖ = pxz[i, k]
                    pxᵢyⱼzₖ = pxyz[i, j, k]
                    f = (pxᵢzₖ / pzₖ) * (pyⱼzₖ / pzₖ)
                    if f != 0.0
                        inner += (pxᵢyⱼzₖ / pzₖ)^q / f^(q-1)
                    end
                end
            end
        end
        cmi += pzₖ * logb(inner)
    end
    return cmi
end
