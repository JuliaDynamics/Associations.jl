export CETsallisFuruichi

"""
    CETsallisFuruichi <: ConditionalEntropy
    CETsallisFuruichi(; base = 2, q = 1.5)

Furuichi (2006)'s discrete Tsallis conditional entropy measure.

## Definition

Furuichi's Tsallis conditional entropy between discrete random variables
``X`` and ``Y`` with finite ranges ``\\mathcal{X}`` and ``\\mathcal{Y}`` is defined as

```math
H_q^T(X | Y) = -\\sum_{x \\in \\mathcal{X}, y \\in \\mathcal{Y}}
p(x, y)^q \\log_q(p(x | y)),
```

when ``q \\neq 1``. For ``q = 1``, ``H_q^T(X | Y)`` reduces to the Shannon conditional
entropy:

```math
H_{q=1}^T(X | Y) = -\\sum_{x \\in \\mathcal{X}, y \\in \\mathcal{Y}} =
p(x, y) \\log(p(x | y))
```
"""
struct CETsallisFuruichi{E} <: ConditionalEntropy
    e::E
    function CETsallisFuruichi(; q = 1.5, base = 2)
        e = MLEntropy(Tsallis(; q, base))
        new{typeof(e)}(e)
    end
end

function estimate(measure::CETsallisFuruichi, c::ContingencyMatrix{T, 2}) where {T}
    e = measure.e.definition
    Nx, Ny = size(c)
    q = e.q
    if q == 1
        return estimate(CEShannon(;base=measure.e.base), pxy)
    end
    py = probabilities(c, dims = 2)
    pxy = probabilities(c)
    ce = 0.0
    qlog = logq0(q)
    for j in 1:Ny
        pyⱼ = py[j]
        for i in 1:Nx
            pxyᵢⱼ = pxy[i, j]
            ce += pxyᵢⱼ^q * qlog(pxyᵢⱼ / pyⱼ)
        end
    end
    ce *= -1.0

    return ce
end


function estimate(measure::CETsallisFuruichi, est::ProbOrDiffEst, x, y)
    throw(ArgumentError("CETsallisFurichi not implemented for $(typeof(est))"))
end
