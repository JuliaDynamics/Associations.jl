export MITsallisMartin

"""
    MITsallisMartin <: MutualInformation
    MITsallisMartin(; base = 2, q = 1.5)

The discrete Tsallis mutual information from Martin et al. (2005)[^Martin2004].

## Usage

- Use with [`independence`](@ref) to perform a formal hypothesis test for pairwise dependence.
- Use with [`mutualinfo`](@ref) to compute the raw mutual information. 

## Description

Martin et al.'s Tsallis mutual information between variables ``X \\in \\mathbb{R}^{d_X}`` and
``Y \\in \\mathbb{R}^{d_Y}`` is defined as

```math
I_{\\text{Martin}}^T(X, Y, q) := H_q^T(X) + H_q^T(Y) - (1 - q) H_q^T(X) H_q^T(Y) - H_q(X, Y),
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

function estimate(measure::MITsallisMartin, est::Contingency{<:ProbabilitiesEstimator}, x...)
    return estimate(measure, contingency_matrix(est.est, x...))
end

function estimate(measure::MITsallisMartin, est::Contingency{<:Nothing}, x...)
    return estimate(measure, contingency_matrix(x...))
end

# This is definition 3 in Martin et al. (2004), but with pᵢ replaced by the joint
# distribution and qᵢ replaced by the product of the marginal distributions.
function estimate(
        measure::MITsallisMartin,
        pxy::ContingencyMatrix{T, 2}) where T
    e = measure.e
    q = measure.e.q
    q != 1 || throw(ArgumentError("MITsallisMartin for q=$(q) not defined with estimator ContingencyMatrix"))
    px = probabilities(pxy, dims = 1)
    py = probabilities(pxy, dims = 2)

    mi = 0.0
    for (i, pxᵢ) in enumerate(px.p)
        for (j, pyⱼ) in enumerate(py.p)
            pxyᵢⱼ = pxy[i, j]
            mi += pxyᵢⱼ^q / (pxᵢ^(q - 1) * pyⱼ^(q - 1))
        end
    end
    f = 1 / (q - 1)
    return f * (1 - mi)
end

function estimate(measure::MITsallisMartin, est::ProbabilitiesEstimator, x, y)
    HX, HY, HXY = marginal_entropies_mi3h(measure, est, x, y)
    q = measure.e.q
    return HX + HY - (1 - q) * HX * HY - HXY
end


function estimate(::MITsallisMartin, est::DifferentialEntropyEstimator, args...)
    throw(ArgumentError("MITsallisMartin not implemented for $(typeof(est))"))
end
