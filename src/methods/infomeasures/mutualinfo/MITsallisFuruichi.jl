export MITsallisFuruichi

"""
    MITsallisFuruichi <: MutualInformation
    MITsallisFuruichi(; base = 2, q = 1.5)

The discrete Tsallis mutual information from Furuichi (2006)[^Furuichi2006], which
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

[^Furuichi2006]:
    Furuichi, S. (2006). Information theoretical properties of Tsallis entropies.
    Journal of Mathematical Physics, 47(2), 023302.

See also: [`mutualinfo`](@ref).
"""
struct MITsallisFuruichi{E <: Tsallis} <: MutualInformation{E}
    e::E
    function MITsallisFuruichi(; q = 1.5, base = 2)
        e = Tsallis(; q, base)
        new{typeof(e)}(e)
    end
end

function estimate(measure::MITsallisFuruichi, est::Contingency{<:ProbabilitiesEstimator}, x...)
    return estimate(measure, contingency_matrix(est.est, x...))
end

function estimate(measure::MITsallisFuruichi, est::Contingency{<:Nothing}, x...)
    return estimate(measure, contingency_matrix(x...))
end

function estimate(
        measure::MITsallisFuruichi,
        pxy::ContingencyMatrix{T, 2}) where T
    e = measure.e
    q = measure.e.q
    px = probabilities(pxy, dims = 1)
    py = probabilities(pxy, dims = 2)

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

function estimate(measure::MITsallisFuruichi, est::ProbabilitiesEstimator, x, y)
    HX, HY, HXY = marginal_entropies_mi3h(measure, est, x, y)
    q = measure.e.q
    return HX + HY - HXY
end


function estimate(::MITsallisFuruichi, est::DifferentialEntropyEstimator, args...)
    throw(ArgumentError("MITsallisFuruichi not implemented for $(typeof(est))"))
end
