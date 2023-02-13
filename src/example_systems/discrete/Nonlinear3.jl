using LabelledArrays: @LArray
using StaticArrays: SVector
using DynamicalSystemsBase: DiscreteDynamicalSystem
using Distributions: Normal
using Random

export Nonlinear3

"""
    Nonlinear3 <: DiscreteDefinition
    Nonlinear3(; xi = rand(3),
        σ₁ = 1.0, σ₂ = 1.0, σ₃ = 1.0,
        a₁ = 3.4, a₂ = 3.4, a₃ = 3.4,
        b₁ = 0.4, b₂ = 0.4, b₃ = 0.4,
        c₁₂ = 0.5, c₂₃ = 0.3, c₁₃ = 0.5,
        rng = Random.default_rng())

A 3d nonlinear system with nonlinear couplings ``x_1 \\to x_2``,
``x_2 \\to x_3`` and ``x_1 \\to x_3``. Modified from Gourévitch et al.
(2006)[Gourévitch2006].

## Equations of motion

```math
\\begin{aligned}
x_1(t+1) &= a_1 x_1 (1-x_1(t))^2  e^{-x_2(t)^2} + 0.4 \\xi_{1}(t) \\\\
x_2(t+1) &= a_1 x_2 (1-x_2(t))^2  e^{-x_2(t)^2} + 0.4 \\xi_{2}(t) + b x_1 x_2 \\\\
x_3(t+1) &= a_3 x_3 (1-x_3(t))^2  e^{-x_3(t)^2} + 0.4 \\xi_{3}(t) + c x_{2}(t) + d x_{1}(t)^2.
\\end{aligned}
```

[Gourévitch2006]:
    Gourévitch, B., Le Bouquin-Jeannès, R., & Faucon, G. (2006). Linear and nonlinear
    causality between signals: methods, examples and neurophysiological
    applications. Biological Cybernetics, 95(4), 349–369.
"""
Base.@kwdef struct Nonlinear3{V, Σx, Σy, Σz, AX, AY, AZ, BX, BY, BZ, C1, C2, C3, RNG} <: DiscreteDefinition
    xi::V = [0.1, 0.2, 0.3]
    σx::Σx = Normal(0, 1.0)
    σy::Σy = Normal(0, 1.0)
    σz::Σz = Normal(0, 1.0)
    ax::AX = 3.4
    ay::AY = 3.4
    az::AZ = 3.4
    bx::BX = 0.4
    by::BY = 0.4
    bz::BZ = 0.4
    c_xy::C1 = 0.5
    c_xz::C2 = 0.3
    c_yz::C3 = 0.5
    rng::RNG = Random.default_rng()
end

function system(definition::Nonlinear3)
    return DiscreteDynamicalSystem(eom_nonlinear3, definition.xi, definition)
end

function eom_nonlinear3(u, p, n)
    x, y, z = u
    (; xi, σx, σy, σz, ax, ay, az, bx, by, bz, c_xy, c_xz, c_yz, rng) = p
    ξ₁ = rand(rng, σx)
    ξ₂ = rand(rng, σy)
    ξ₃ = rand(rng, σz)
    dx = ax*x*(1-x)^2 * exp(-x^2) + bx*ξ₁
    dy = ay*y*(1-y)^2 * exp(-y^2) + by*ξ₂ + c_xy*x*y
    dz = az*z*(1-z)^2 * exp(-z^2) + bz*ξ₃ + c_yz*y + c_xz*x^2
    return SVector{3}(dx, dy, dz)
end
