export repressilator

using LabelledArrays: @LArray
using StaticArrays: SVector
using DynamicalSystemsBase: trajectory
using DynamicalSystemsBase: ContinuousDynamicalSystem

export Repressilator6

"""
    Repressilator6 <: ContinuousDefinition
    Repressilator6(; xi = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6], α = 10.0, α₀ = 0.0, β = 100.0,
        n = 2) → ContinuousDynamicalSystem

A six-dimensional repressilator (or repression-driven oscillator) from Elowitz & Leibler
(2000)[^Elowitz2000].

Used in Sun & Bollt (2014)[^Sun2014] to study the performance of the causation entropy
algorithm.

## Description

```math
\\begin{align*}
\\dfrac{dm_1}{dt} &= -m1 + \\dfrac{\\alpha}{1 + p_3^n} + \\alpha_0 \\\\
\\dfrac{dm_2}{dt} &= -m2 + \\dfrac{\\alpha}{1 + p_1^n} + \\alpha_0 \\\\
\\dfrac{dm_3}{dt} &= -m3 + \\dfrac{\\alpha}{1 + p_2^n} + \\alpha_0 \\\\
\\dfrac{dp_1}{dt} &= -\\beta(p_1 - m_1) \\\\
\\dfrac{dp_2}{dt} &= -\\beta(p_2 - m_2) \\\\
\\dfrac{dp_3}{dt} &= -\\beta(p_3 - m_3) \\\\
\\end{align*}
```

[^Elowitz2000]: Elowitz, M. B., & Leibler, S. (2000). A synthetic oscillatory network of
    transcriptional regulators. Nature, 403(6767), 335-338.
[^Sun2014]: Sun, J., Cafaro, C., & Bollt, E. M. (2014). Identifying the coupling structure
    in complex systems through the optimal causation entropy principle. Entropy, 16(6),
    3416-3433.
"""
Base.@kwdef struct Repressilator6{V, A, A0, B, N} <: ContinuousDefinition
    xi::V = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    α::A = 10.0
    α₀::A0 = 0.0
    β::B = 100.0
    n::N = 2
end

function system(definition::Repressilator6)
    return ContinuousDynamicalSystem(eom_repressilator6, definition.xi, definition)
end

@inline @inbounds function eom_repressilator6(u, p, t)
    (; xi, α, α₀, n, β) = p
    # pᵢ := concentration of protein repressor i
    # mᵢ := concentration of mRNA associated with pᵢ
    m₁, m₂, m₃, p₁, p₂, p₃ = u

    ṁ₁ = -m₁ + α/(1 + p₃^n) + α₀
    ṁ₂ = -m₂ + α/(1 + p₁^n) + α₀
    ṁ₃ = -m₃ + α/(1 + p₂^n) + α₀
    ṗ₁ = -β*(p₁ - m₁)
    ṗ₂ = -β*(p₂ - m₂)
    ṗ₃ = -β*(p₃ - m₃)

    return SVector{6}(ṁ₁, ṁ₂, ṁ₃, ṗ₁, ṗ₂, ṗ₃)
end
