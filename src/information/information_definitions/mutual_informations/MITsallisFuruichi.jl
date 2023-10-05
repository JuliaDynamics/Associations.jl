using ComplexityMeasures: Tsallis

export MITsallisFuruichi
"""
    MITsallisFuruichi <: BivariateInformationMeasure
    MITsallisFuruichi(; base = 2, q = 1.5)

The discrete Tsallis mutual information from Furuichi (2006)[Furuichi2006](@cite), which
in that paper is called the *mutual entropy*.

## Usage

- Use with [`independence`](@ref) to perform a formal hypothesis test for pairwise dependence.
- Use with [`mutualinfo`](@ref) to compute the raw mutual information.


## Description

Furuichi's Tsallis mutual entropy between variables ``X \\in \\mathbb{R}^{d_X}`` and
``Y \\in \\mathbb{R}^{d_Y}`` is defined as

```math
I_q^T(X; Y) = H_q^T(X) - H_q^T(X | Y) = H_q^T(X) + H_q^T(Y) - H_q^T(X, Y),
```

where ``H^T(\\cdot)`` and ``H^T(\\cdot, \\cdot)`` are the marginal and joint Tsallis
entropies, and `q` is the [`Tsallis`](@ref)-parameter.
```

See also: [`mutualinfo`](@ref).
"""
struct MITsallisFuruichi{E <: Tsallis} <: MutualInformation
    e::E
    function MITsallisFuruichi(; q = 1.5, base = 2)
        e = Tsallis(; q, base)
        new{typeof(e)}(e)
    end
end

function information(definition::MITsallisFuruichi, pxy::Probabilities{T, 2}) where T
    e = definition.e
    q = definition.e.q
    px = marginal(pxy, dims = 1)
    py = marginal(pxy, dims = 2)

    mi = 0.0
    for i in eachindex(px.p)
        for j in eachindex(py.p)
            pxyᵢⱼ = pxy[i, j]
            mi += pxyᵢⱼ^q / (px[i]^(q - 1) * py[j]^(q - 1))
        end
    end
    mi = (1 / (q - 1) * (1 - mi) / (1-q))
    return _convert_logunit(mi, ℯ, e.base)
end