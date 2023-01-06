export MITsallisFuruichi

"""
    MITsallisFuruichi <: MutualInformation
    MITsallisFuruichi(; base = 2, q = 1.5)


The mutual *entropy* formulation of Tsallis "mutual information"
(Furuichi, 2006)[^Furuichi2006], which is the direct analogue of the 3-entropies
decomposition of Shannon mutual information ([`MIShannon`](@ref)=.

## Description

Furuichi's Tsallis mutual entropy between variables ``X \\in \\mathbb{R}^{d_X}`` and
``Y \\in \\mathbb{R}^{d_Y}`` is defined as

```math
I_q^T(X; Y) = H_q^T(X) - H_q^T(X | Y) = H_q^T(X) + H_q^T(Y) - H_q^T(X, Y),
```

where ``H^T(\\cdot)`` and ``H^T(\\cdot, \\cdot)`` are the marginal and joint Tsallis
entropies, and `q` is the [`Tsallis`](@ref)-parameter. Equivalently,
Vilasini & Colbeck (2019) give the direct definition (compatible with
[`ContingencyMatrix`](@ref))

```math
I(X, Y; Z)^R_q =
\\dfrac{1}{q-1}
\\log \\left(
    \\sum_{x \\in X}\\sum_{y \\in Y}
    \\dfrac{p(x, y)^q}{\\left( p(x)\\cdot p(y) \\right)^{q-1}}
\\right)
```

[^Furuichi2006]:
    Furuichi, S. (2006). Information theoretical properties of Tsallis entropies.
    Journal of Mathematical Physics, 47(2), 023302.
[^Vilasini2019]: Vilasini, V., & Colbeck, R. (2019). Analyzing causal structures using
    Tsallis entropies. Physical Review A, 100(6), 062108.

See also: [`mutualinfo`](@ref).
"""
struct MITsallisFuruichi{E <: Tsallis} <: MutualInformation{E}
    e::E
    function MITsallisFuruichi(; q = 1.5, base = 2)
        e = Tsallis(; q, base)
        new{typeof(e)}(e)
    end
end

function estimate(
        measure::MITsallisFuruichi,
        pxy::ContingencyMatrix{T, 2}) where T
    e = measure.e
    q = measure.e.q
    px = probabilities(pxy, 1)
    py = probabilities(pxy, 2)

    mi = 0.0
    for i in eachindex(px.p)
        for j in eachindex(py.p)
            pxyᵢⱼ = pxy[i, j]
            mi += pxyᵢⱼ^q / (px[i]^(q - 1) * py[j]^(q - 1))
        end
    end
    return (1 / (q - 1) * (1 - mi) / (1-q)) / log(e.base, ℯ)
end

function estimate(measure::MITsallisFuruichi, est::ProbOrDiffEst, x, y)
    e = measure.e
    X = Dataset(x); Y = Dataset(y); XY = Dataset(X, Y)
    HX = entropy(e, est, X)
    HY = entropy(e, est, Y)
    HXY = entropy(e, est, XY)
    return HX + HY - HXY
end

function estimate(measure::MITsallisFuruichi, est::WellDefinedMIShannonProbEsts{m, D},
        x, y) where { m, D}
    e = measure.e
    pX, pY, pXY = marginal_probabilities(measure, est, x, y)
    e = measure.e
    HX = entropy(e, pX)
    HY = entropy(e, pY)
    HXY = entropy(e, pXY)
    return HX + HY - HXY
end
