using Random
using Distributions: Normal
using DynamicalSystemsBase: ContinuousDynamicalSystem
using StaticArrays: SVector

export ChuaCircuitsBidir6

"""
    ChuaCircuitsBidir6 <: ContinuousDefinition
    ChuaCircuitsBidir6(;u₀ = [0.1, 0.1, 0.2, 0.15, 0.15, 0.22],
        α₁ = 7.0, α₂ = 7.0, β₁ = 14.286, β₂ = 14.286,
        F₁ = 1.5, F₂ = 1.5, ω₁ = 3.0, ω₂ = 3.0,
        n1 = Normal(0, 0.1),
        n2 = Normal(0, 0.1),
        c12 = 0.1, c21 = 0.1, m₀ = -1/7, m₁ = 2/7)

Initialize a bidirectionally coupled system consisting of two driven Chua
circuits, X₁ and X₂ (Murali &  Lakshmanan, 1993)[^Murali1993].

## Description

The subsystems are mutually coupled by a linear resistor, where `ϵ12` controls the
influence of X₁ on X₂, and `c21` controls the influence of X₂ on X₁. The parameters for
the subsystems are set equal to each other, as in the original paper, but can here
be tuned individually for each subsystem.

```math
\\begin{align*}
\\dfrac{dx_1}{dt} &= \\alpha_1(y_1, h(x_1)) - \\alpha_1 \\epsilon(x_1 - x_2) \\\\
\\dfrac{dy_1}{dt} &= x_1 - y_1 + z_1 \\\\
\\dfrac{dz_1}{dt} &= -\\beta_1 y_1 + F_1 sin(\\omega_1 t) + \\epsilon_1 \\\\
\\dfrac{dx_2}{dt} &= \\alpha_2 (y_2, h(x_2)) - \\alpha_2 c_{12}(x_1 - x_2) \\\\
\\dfrac{dy_2}{dt} &= x_2 - y_2 + z_2 \\\\
\\dfrac{dz_2}{dt} &= -\\beta_2 y_2 + F_2 sin(\\omega_2 t) + \\epsilon_2 \\\\,
\\end{align*}
```

where ``h(x) = M_1x + 0.5(M_0 - M_1)(|x+1| - |x - 1|)`` and ``\\epsilon_1, \\epsilon_2``
are noise terms that at each integration step is drawn independently from the normal
distributions `n1` and `n2`, respectively.

[^Murali1993]:
    Murali, K., and M. Lakshmanan. "Chaotic dynamics of the driven Chua's
    circuit." IEEE Transactions on Circuits and Systems I Fundamental
    Theory and Applications 40.11 (1993): 836-840.
"""
Base.@kwdef struct ChuaCircuitsBidir6{V,A1,A2,B1,B2,F1,F2,W1,W2,S1,S2,E1,E2,M0,M1,R}<: ContinuousDefinition
    xi::V = [0.1, 0.1, 0.2, 0.15, 0.15, 0.22]
    α₁::A1 = 7.0
    α₂::A2 = 7.0
    β₁::B1 = 14.286
    β₂::B2 = 14.286
    F₁::F1 = 1.5
    F₂::F2 = 1.5
    ω₁::W1 = 3.0
    ω₂::W2 = 3.0
    n1::S1 = Normal(0, 0.1)
    n2::S2 = Normal(0, 0.1)
    c12::E1 = 0.1
    c21::E2 = 0.1
    m₀::M0 = -1/7
    m₁::M1 = 2/7
    rng::R = Random.default_rng()
end

function system(definition::ChuaCircuitsBidir6)
    return ContinuousDynamicalSystem(eom_chuabidir6, definition.xi, definition)
end

function eom_chuabidir6(u, p::ChuaCircuitsBidir6, t)
    (; xi, α₁, α₂, β₁, β₂, F₁, F₂, ω₁, ω₂, n1, n2, c12, c21, m₀, m₁, rng) = p
    x₁, y₁, z₁,  x₂, y₂, z₂ = u

    ξ1 = rand(rng, n1)
    ξ2 = rand(rng, n2)
    hx₁ = m₁*x₁ + 0.5*(m₀ - m₁)*(abs(x₁+1) - abs(x₁-1))
    hx₂ = m₁*x₂ + 0.5*(m₀ - m₁)*(abs(x₂+1) - abs(x₂-1))

    dx₁ = α₁*(y₁-hx₁) - α₁*c21*(x₁ - x₂)
    dy₁ = x₁-y₁+z₁
    dz₁ = -β₁*y₁ + F₁*sin(ω₁*t) + ξ1

    dx₂ = α₂*(y₂-hx₂) - α₂*c12*(x₁ - x₂)
    dy₂ = x₂-y₂+z₂
    dz₂ = -β₂*y₂ + F₂*sin(ω₂*t) + ξ2
    SVector{6}(dx₁, dy₁, dz₁, dx₂, dy₂, dz₂)
end
