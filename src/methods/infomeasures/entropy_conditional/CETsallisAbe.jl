export CETsallisAbe

"""
    CETsallisAbe <: ConditionalEntropy
    CETsallisAbe(; base = 2, q = 1.5)

Abe & Rajagopal (2001)'s discrete Tsallis conditional entropy measure.

## Definitions

Abe & Rajagopal's Tsallis conditional entropy between discrete random variables
``X`` and ``Y`` with finite ranges ``\\mathcal{X}`` and ``\\mathcal{Y}`` is defined as

```math
H_q^{T_A}(X | Y) = \\dfrac{1}{1-q}
\\dfrac{\\sum_{x \\in \\mathcal{X}, y \\in \\mathcal{Y}} = p(x, y)^q \\log_q(p(x | y))}{\\sum_{y \\in \\mathcal{Y}} p(y)^q},
```

[^Abe2001]:
    Abe, S., & Rajagopal, A. K. (2001). Nonadditive conditional entropy and its
    significance for local realism. Physica A: Statistical Mechanics and its Applications, 289(1-2), 157-164.
"""
struct CETsallisAbe{E} <: ConditionalEntropy
    e::E
    function CETsallisAbe(; q = 1.5, base = 2)
        e = Tsallis(; q, base)
        new{typeof(e)}(e)
    end
end

function estimate(measure::CETsallisAbe, pxy::ContingencyMatrix{T, 2}) where {T}
    e = measure.e
    Nx, Ny = size(pxy)
    base, q = e.base, e.q

    py = probabilities(pxy, 2)
    ce = 0.0
    qlog = logq0(q)
    for j in 1:Ny
        pyⱼ = py[j]
        for i in 1:Nx
            pxyᵢⱼ = pxy[i, j]
            ce += pxyᵢⱼ^q * qlog(pxyᵢⱼ / pyⱼ)
        end
    end
    ce *= -1
    ce /= sum(py .^ q)

    if q == 1
        return (1/(1 - q))* (ce - 1) / log(base, ℯ)
    else
        return -(1/(1 - q))* (ce - 1)
    end
end
