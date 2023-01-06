export CEShannon

"""
    CEShannon <: ConditionalEntropy
    CEShannon(; base = 2,)

The [`Shannon`](@ref) conditional entropy measure.

## Definition

The conditional entropy between discrete random variables
``X`` and ``Y`` with finite ranges ``\\mathcal{X}`` and ``\\mathcal{Y}`` is defined as

```math
-\\sum_{x \\in \\mathcal{X}, y \\in \\mathcal{Y}} = p(x, y) \\log(p(x | y)).
```
"""
struct CEShannon{E} <: ConditionalEntropy
    e::E
    function CEShannon(; base = 2)
        e = Shannon(; base)
        new{typeof(e)}(e)
    end
end

function estimate(measure::CEShannon, pxy::ContingencyMatrix{T, 2}) where {T}
    e = measure.e
    Nx, Ny = size(pxy)
    py = probabilities(pxy, 2)
    ce = 0.0
    log0 = log_with_base(e.base)
    for j in 1:Ny
        pyⱼ = py[j]
        for i in 1:Nx
            pxyᵢⱼ = pxy[i, j]
            ce += pxyᵢⱼ * log0(pxyᵢⱼ / pyⱼ)
        end
    end
    return -ce
end
