using ComplexityMeasures: Renyi

export CMIRenyiPoczos

"""
    CMIRenyiPoczos <: ConditionalMutualInformation
    CMIRenyiPoczos(; base = 2, q = 1.5)

The differential Rényi conditional mutual information ``I_q^{R_{P}}(X; Y | Z)``
defined in [Poczos2012](@citet).

## Usage

- Use with [`association`](@ref) to compute the raw Rényi-Poczos conditional mutual information
    using of of the estimators listed below.
- Use with [`independence`](@ref) to perform a formal hypothesis test for pairwise conditional 
    independence using the Rényi-Poczos conditional mutual information.

## Compatible estimators

- [`PoczosSchneiderCMI`](@ref)

## Definition

```math
\\begin{align*}
I_q^{R_{P}}(X; Y | Z) &= \\dfrac{1}{q-1}
\\int \\int \\int \\dfrac{p_Z(z) p_{X, Y | Z}^q}{( p_{X|Z}(x|z) p_{Y|Z}(y|z) )^{q-1}} \\\\
&= \\mathbb{E}_{X, Y, Z} \\sim p_{X, Y, Z}
\\left[ \\dfrac{p_{X, Z}^{1-q}(X, Z) p_{Y, Z}^{1-q}(Y, Z) }{p_{X, Y, Z}^{1-q}(X, Y, Z) p_Z^{1-q}(Z)} \\right]
\\end{align*}
```

## Estimation

- [Example 1](@ref CMIRenyiPoczos_PoczosSchneiderCMI): Dedicated [`PoczosSchneiderCMI`](@ref) estimator.
"""
Base.@kwdef struct CMIRenyiPoczos{B, Q} <: ConditionalMutualInformation
    base::B = 2
    q::Q = 1.5
end
