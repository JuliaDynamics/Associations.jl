export MIDefinitionTsallisH3Martin

"""
    MIDefinitionTsallisH3Martin <: MutualInformationDefinition

A directive to compute discrete Tsallis mutual information
(Martin et al., 2005)[^Martin2004].

Used in combination with [`MITsallis`](@ref) and [`mutualinfo`](@ref).

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
"""
struct MIDefinitionTsallisH3Martin <: MIH3 end
