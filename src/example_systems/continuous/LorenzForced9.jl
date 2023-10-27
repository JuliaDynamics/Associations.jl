using LabelledArrays: @LArray
using StaticArrays: SVector
using DynamicalSystemsBase: trajectory
using DynamicalSystemsBase: ContinuousDynamicalSystem
using Distributions: Uniform

export LorenzForced9

"""
    LorenzForced9{V} <: ContinuousDefinition
    LorenzForced9(; xi = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
        c_xy = 0.1, c_yx = 0.1,
        c_zx = 0.05, c_zy = 0.05,
        a₁ = 10, a₂ = 28, a₃ = 8/3,
        b₁ = 10, b₂ = 28, b₃ = 8/3,
        c₁ = 10, c₂ = 28, c₃ = 8/3)

A system consisting of two bidirectionally coupled 3D Lorenz
systems forced by an external 3D Lorenz system (Amigó & Hirata, 2018).

## Description

The dynamics is generated by the following vector field

```math
\\begin{align*}
\\dot{x_1} &= - a_1 (x_1 - x_2) + c_{yx}(y_1 - x_1) + c_{zx}(z_1 - x_1) \\\\
\\dot{x_2} &= - x_1 x_3 + a_2 x_1 - x_2 \\\\
\\dot{x_3} &= x_1 x_2 - a_3 x_3 \\\\
\\dot{y_1} &= -b_1 (y_1 - y_2) + c_{xy} (x_1 - y_1) + c_{zy}(z_1 - y_1) \\\\
\\dot{y_2} &= - y_1 y_3 + b_2 y_1 - y_2 \\\\
\\dot{y_3} &= y_1 y_2 - b_3 y_3 \\\\
\\dot{z_1} &= - c_1 (z_1 - z_2) \\\\
\\dot{z_2} &= - z_1 z_3 + c_2 z_1 - z_2 \\\\
\\dot{z_3} &= z_1 z_2 - c_3 z_3
\\end{align*}
```

[^Amigó2018]:
    Amigó, José M., and Yoshito Hirata. "Detecting directional couplings from
    multivariate flows by the joint distance distribution." Chaos: An
    Interdisciplinary Journal of Nonlinear Science 28.7 (2018): 075302.
"""
Base.@kwdef struct LorenzForced9{V,CXY,CYX,CZX,CZY,A1,A2,A3,B1,B2,B3,C1,C2,C3} <: ContinuousDefinition
    xi::V = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    c_xy::CXY = 1.0
    c_yx::CYX = 1.0
    c_zx::CZX = 1.0
    c_zy::CZY = 1.0 # beyond c = 2, systems synchronize
    a₁::A1 = 10
    a₂::A2 = 28
    a₃::A3 = 8/3
    b₁::B1 = 10
    b₂::B2 = 28
    b₃::B3 = 8/3
    c₁::C1 = 10
    c₂::C2 = 28
    c₃::C3 = 8/3
end

function system(definition::LorenzForced9)
    return ContinuousDynamicalSystem(eom_lorenzforced9, definition.xi, definition)
end

@inline @inbounds function eom_lorenzforced9(u, p::LorenzForced9, t)
    (; xi, c_xy, c_yx, c_zx, c_zy, a₁, a₂, a₃, b₁, b₂, b₃, a₃, c₁, c₂, c₃) = p
    x₁, x₂, x₃, y₁, y₂, y₃, z₁, z₂, z₃ = u

    dx₁ = -a₁*(x₁ - x₂) + c_yx*(y₁ - x₁) + c_zx*(z₁ - x₁)
    dx₂ = -x₁*x₃ + a₂*x₁ - x₂
    dx₃ = x₁*x₂ - a₃*x₃

    dy₁ = -b₁*(y₁ - y₂) + c_xy*(x₁ - y₁) + c_zy*(z₁ - y₁)
    dy₂ = -y₁*y₃ + b₂*y₁ - y₂
    dy₃ = y₁*y₂ - b₃*y₃

    dz₁ = -c₁*(z₁ - z₂)
    dz₂ = -z₁*z₃ + c₂*z₁ - z₂
    dz₃ = z₁*z₂ - c₃*z₃

    return SVector{9}(dx₁, dx₂, dx₃, dy₁, dy₁, dy₃, dz₁, dz₂, dz₃)
end
