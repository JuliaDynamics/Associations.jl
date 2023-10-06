export CETsallisFuruichi

"""
    CETsallisFuruichi <: ConditionalEntropy
    CETsallisFuruichi(; base = 2, q = 1.5)

Furuichi (2006)'s discrete Tsallis conditional entropy definition.

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

If any of the entries of the marginal distribution for `Y` are zero, then the 
measure is undefined and `NaN` is returned.
"""
struct CETsallisFuruichi{E} <: ConditionalEntropy
    e::E
    function CETsallisFuruichi(; q = 1.5, base = 2)
        e = Tsallis(; q, base)
        new{typeof(e)}(e)
    end
end

function information(definition::CETsallisFuruichi, pxy::Probabilities{T, 2}) where {T}
    e = definition.e
    Nx, Ny = size(pxy)
    q = e.q
    if q == 1
        return information(CEShannon(; base = definition.e.base), pxy)
    end
    py = marginal(pxy, dims = 2)
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


function logq0(q)
    if q == 1.0
        return x -> zero(x)
    else
        return x -> (x^(1 - q) - 1)/(1 - q)
    end
end

