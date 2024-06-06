using ComplexityMeasures: Renyi

export CMIRenyiPoczos

"""
    CMIRenyiPoczos <: ConditionalMutualInformation
    CMIRenyiPoczos(; base = 2, q = 1.5)

The differential Rényi conditional mutual information ``I_q^{R_{P}}(X; Y | Z)``
defined in (Póczos & Schneider, 2012)[^Póczos2012].

## Usage

- Use with [`independence`](@ref) to perform a formal hypothesis test for pairwise dependence.
- Use with [`condmutualinfo`](@ref) to compute the raw conditional mutual information.

## Definition

```math
\\begin{align*}
I_q^{R_{P}}(X; Y | Z) &= \\dfrac{1}{q-1}
\\int \\int \\int \\dfrac{p_Z(z) p_{X, Y | Z}^q}{( p_{X|Z}(x|z) p_{Y|Z}(y|z) )^{q-1}} \\\\
&= \\mathbb{E}_{X, Y, Z} \\sim p_{X, Y, Z}
\\left[ \\dfrac{p_{X, Z}^{1-q}(X, Z) p_{Y, Z}^{1-q}(Y, Z) }{p_{X, Y, Z}^{1-q}(X, Y, Z) p_Z^{1-q}(Z)} \\right]
\\end{align*}
```

[^Póczos2012]:
    Póczos, B., & Schneider, J. (2012, March). Nonparametric estimation of conditional
    information and divergences. In Artificial Intelligence and Statistics (pp. 914-923).
    PMLR.
"""
Base.@kwdef struct CMIRenyiPoczos{B, Q} <: ConditionalMutualInformation
    base::B = 2
    q::Q = 1.5
end
