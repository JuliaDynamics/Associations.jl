export CMITsallisFuruichi

"""
    CMIDefinitionTsallisH4 <: MutualInformationDefinition
    CMIDefinitionTsallisH4Furuichi()

The 4-entropies formulation of Tsallis conditional mutual information (CMI)
(Vilasini & Colbeck, 2019), which is a conditional version of the Tsallis mutual *entropy*
formulation of Tsallis mutual information (Furuichi, 2006)[^Furuichi2006].

Used in combination with [`CMITsallis`](@ref) and [`condmutualinfo`](@ref).

## Description

Vilanisi & Colbeck's Tsallis CMI for random variable  ``X \\in \\mathbb{R}^{d_X}``,
``Y \\in \\mathbb{R}^{d_Y}`` and ``Z \\in \\mathbb{R}^{d_Z}`` is defined as

```math
I_q^T(X; Y | Z) = H_q^T(XZ) + H_q^T(YZ) - H_q^T(X, Y, Z) - H_q^T(Z)
```

where ``H^T(\\cdot)`` and ``H^T(\\cdot, \\cdot)`` are the marginal and joint Tsallis
entropies, and `q` is the [`Tsallis`](@ref)-parameter.

[^Furuichi2006]:
    Furuichi, S. (2006). Information theoretical properties of Tsallis entropies.
    Journal of Mathematical Physics, 47(2), 023302.

[^Vilasini2019]: Vilasini, V., & Colbeck, R. (2019). Analyzing causal structures using
    Tsallis entropies. Physical Review A, 100(6), 062108.

See also: [`mutualinfo`](@ref).
"""
struct CMITsallisFuruichi <: MIH3 end
