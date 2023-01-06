export MITsallisMartin

"""
    MITsallisMartin <: MutualInformation
    MITsallisMartin(; base = 2, q = 1.5)

The discrete Tsallis mutual information from Martin et al. (2005)[^Martin2004].

## Description

Martin et al.'s Tsallis mutual information between variables ``X \\in \\mathbb{R}^{d_X}`` and
``Y \\in \\mathbb{R}^{d_Y}`` is defined as

```math
I_{\\text{Martin}^T(X, Y, q) := H_q^T(X) + H_q^T(Y) - (1 - q)H_q^T(X)H_q^T(Y) - H_q(X, Y),
```

where ``H^S(\\cdot)`` and ``H^S(\\cdot, \\cdot)`` are the marginal and joint Shannon
entropies, and `q` is the [`Tsallis`](@ref)-parameter.

[^Martin2004]:
    Martin, S., Morison, G., Nailon, W., & Durrani, T. (2004). Fast and accurate image
    registration using Tsallis entropy and simultaneous perturbation stochastic
    approximation. Electronics Letters, 40(10), 1.

See also: [`mutualinfo`](@ref).
"""
struct MITsallisMartin{E <: Tsallis} <: MutualInformation{E}
    e::E
    function MITsallisMartin(; q = 1.5, base = 2)
        e = Tsallis(; q, base)
        new{typeof(e)}(e)
    end
end

# function estimate(
#         measure::MITsallisMartin,
#         pxy::ContingencyMatrix{T, 2}) where T
#     e = measure.e
#     q = measure.e.q
#     px = probabilities(pxy, 1)
#     py = probabilities(pxy, 2)

#     mi = 0.0
#     for i in eachindex(px.p)
#         for j in eachindex(py.p)
#             pxyᵢⱼ = pxy[i, j]
#             mi += pxyᵢⱼ^q / (px[i]^(q - 1) * py[j]^(q - 1))
#         end
#     end
#     return (1 / (q - 1) * (1 - mi) / (1-q)) / log(e.base, ℯ)
# end

function estimate(measure::MITsallisMartin, est::ProbOrDiffEst, x, y)
    HX, HY, HXY = marginal_entropies_mi3h(measure, est, x, y)
    q = measure.e.q
    return HX + HY - (1 - q) * HX * HY - HXY
end
