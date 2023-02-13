using DynamicalSystemsBase: ContinuousDynamicalSystem
using StaticArrays: SVector

export RosslerBidir6

"""
    RosslerBidir6 <: ContinuousDefinition
    RosslerBidir6(; xi = [0.1, 0.1, 0.2, 0.3, 0.3, 0.4],
        a = 0.1, b = 0.1, c = 14.0, ϵ₁ = 0.0, ϵ₂ = 0.0,
        ω₁ = 1 + 0.015, ω₂ = 1 - 0.015)

A bidirectionally coupled 6D Rossler system from Krakovská et al. (2018)[^Krakovská2018].

## Description

The system consists of two separate subsystems, each being a 3D Rossler
attractor. The subsystems are bidirectionally coupled, influencing each other
through variables ``x_1`` and ``x_2`.

```math
\\begin{align*}
\\dfrac{dx_1}{dt} = \\omega_1 (-y_1) - z_1 + c_{21}*(x_1 - x_2) \\\\
\\dfrac{dy_1}{dt}  = \\omega_1 x_1 + a y_1 \\\\
\\dfrac{dz_1}{dt}  = b + z_1 (x_1 - c) \\\\
\\dfrac{dx_2}{dt}  = \\omega_2 (-y_2) - z_2 + c_{12} (x_2 - x_1) \\\\
\\dfrac{dy_2}{dt} = \\omega_2*x_2 + a*y_2 \\\\
\\dfrac{dz_2}{dt} = b + z_2 (x_2 - c) \\\\
\\end{align*}
```

with ``c_{12} \\geq 0`` and ``c \\geq 0``.

[^Krakovská2018]:
    Krakovská, A., Jakubík, J., Chvosteková, M., Coufal, D., Jajcay, N., & Paluš, M.
    (2018). Comparison of six methods for the detection of causality in a bivariate time
    series. Physical Review E, 97(4), 042207.
"""
Base.@kwdef struct RosslerBidir6{V, A, B, C, E1, E2, Ω1, Ω2} <: ContinuousDefinition
    xi::V = [0.1, 0.1, 0.2, 0.3, 0.3, 0.4]
    a::A = 0.1
    b::B = 0.1
    c::C = 14.0
    c12::E1 = 0.0
    c21::E2 = 0.0
    ω₁::Ω1 = 1 + 0.015
    ω₂::Ω2 = 1 - 0.015
end

function system(definition::RosslerBidir6)
    return ContinuousDynamicalSystem(eom_rosslerrosslerbidir6, definition.xi, definition)
end

function eom_rosslerrosslerbidir6(u, p, t)
    (; xi, a, b, c, c12, c21, ω₁, ω₂) = p
    x₁, y₁, z₁, x₂, y₂, z₂ = u

    # First Rössler system
    dx₁ = ω₁*(-y₁) - z₁ + c21*(x₁ - x₂)
    dy₁ = ω₁*x₁ + a*y₁
    dz₁ = b + z₁*(x₁ - c)
    # Second Rössler system
    dx₂ = ω₂*(-y₂) - z₂ + c12*(x₂ - x₁)
    dy₂ = ω₂*x₂ + a*y₂
    dz₂ = b + z₂*(x₂ - c)
    return SVector{6}(dx₁, dy₁, dz₁, dx₂, dy₂, dz₂)
end
