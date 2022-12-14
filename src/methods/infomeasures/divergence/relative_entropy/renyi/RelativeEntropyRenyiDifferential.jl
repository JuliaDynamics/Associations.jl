
export RelativeEntropyRenyiDifferential
export RenyiDivergenceDifferential

"""
    RenyiDivergenceDifferential  <: DivergenceDefinition
    RenyiDivergenceDifferential()

An instruction to compute the Rényi divergence (relative entropy) according to the
original definition in Rényi (1961)[^Rényi1961]'s seminal paper, here stated in terms of
notation from Van Erven et al. (2014)[[^VanErven2014]].

## Description

For a discrete sample space ``\\Omega`` and probability mass functions
``p(x) : \\Omega \\to [0, 1]`` and ``q(x) : \\Omega \\to [0, 1]``,
the Rényi relative entropy (divergence) is, for `q != 1`, given by

```math
D_q(P || Q) = \\dfrac{1}{q - 1} \\log \\sum_{i = 1}^n p_i^q q_i^{1-\\alpha}.
```

[^Rényi1961]:
    Rényi, A. (1961, June). On measures of entropy and information. In Proceedings of the
    fourth Berkeley symposium on mathematical statistics and probability (Vol. 1, No.
    547-561).

[^VanErven2014]:
    Van Erven, T., & Harremos, P. (2014). Rényi divergence and Kullback-Leibler divergence.
    IEEE Transactions on Information Theory, 60(7), 3797-3820.
"""
struct RenyiDivergenceDifferential <: DivergenceDefinition end

"""
    RelativeEntropyRenyiDifferential <: MutualInformation
    RelativeEntropyRenyiDifferential(; base = 2, q = 1.5,
        definition = RenyiDivergenceDifferential())

`RelativeEntropyRenyiDifferential` is a directive to compute the discrete Rényi relative entropy
(or divergence) to the given `base` between random variables
``X \\in \\mathbb{R}^{d_X}`` and ``Y \\in \\mathbb{R}^{d_Y}`` using the formula.

## Supported definitions

- [`RenyiDivergenceDifferential`](@ref).

See also: [`divergence`](@ref).
"""
struct RelativeEntropyRenyiDifferential{D <: DivergenceDefinition, E <: Entropy} <: Divergence
    e::E
    definition::D
    function RelativeEntropyRenyiDifferential(; base = 2, q = 1.5,
            definition::D = RenyiDivergenceDifferential()) where {D}
            e = Renyi(; base, q)
        new{D, typeof(e)}(e, definition)
    end
end


function estimate(e::RelativeEntropyRenyiDifferential{<:RenyiDivergenceDifferential},
        est::ProbabilitiesEstimator, x, y)
    q = e.q
    p = probabilities(est.est, x)
    q = probabilities(est.est, y)
    re = 1 / (q - 1) * log(sum(pᵢ^q * qᵢ^(1 - q) for (pᵢ, qᵢ) in zip(p, q)))
    return re / log(e.base, ℯ)
end


# This function contain analytical expressions for various relative entropies.
# Renyi divergence expressions are from Gil, M. (2011). On Rényi divergence measures for continuous alphabet sources. PhD Thesis.
using Distributions: Beta
using SpecialFunctions: gamma
_beta(x, y) = gamma(x)*gamma(y) / gamma(x + y)

function divergence(e::RelativeEntropyRenyiDifferential, x::Beta, y::Beta)
    q = e.q # Our q is their α
    αx, αy, βx, βy = x.α, y.α, x.β, y.β
    a = q*αx + (1 - q)*αy
    b = q*βx + (1 - q)*βy
    @assert a >= 0 && b >= 0
    re = log(_beta(αy, βy) / _beta(αx, βx)) +
        (1 / (q - 1)) * log(_beta(a, b) / _beta(αx, βx))
    return re / log(e.base, ℯ)
end


# Analytical KL divergence for 1D normals:
# https://ieeexplore.ieee.org/abstract/document/6832827?casa_token=ZhfFH5_G6XgAAAAA:RzQMg0Zjn-CwtOWw4N-jeum3bWzP7tRioSSFAb76fZX58JmXDBW7mSqjbxvr73NDa9fplUSIGw
