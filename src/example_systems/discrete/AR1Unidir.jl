using LabelledArrays: @LArray
using StaticArrays: SVector
using DynamicalSystemsBase: DiscreteDynamicalSystem
using Distributions: Normal
using Random

export AR1Unidir

"""
    AR1Unidir <: DiscreteDefinition
    AR1Unidir(; ui = [0.2, 0.3], a₁ = 0.90693, b₁ = 0.40693, c_xy = 0.5,
        nx = Normal(0, 0.40662), ny = Normal(0, 0.40662),
        rng::R = Random.default_rng())

A bivariate, order one autoregressive model, where ``x \\to y`` (Paluš et al,
2018)[^Paluš2018].

## Equations of motion

```math
\\begin{aligned}
x(t+1) &= a_1 x(t) + \\xi_{1} \\\\
y(t+1) &= b_1 y(t) - c_{xy} x + \\xi_{2},
\\end{aligned}
```

where ``\\xi_{1}`` and ``\\xi_{2}`` are drawn from normal distributions `nx` and `ny`
at each iteration.

[^Paluš2018]:
    Paluš, M., Krakovská, A., Jakubík, J., & Chvosteková, M. (2018). Causality,
    dynamical systems and the arrow of time. Chaos: An Interdisciplinary Journal of
    Nonlinear Science, 28(7), 075307. http://doi.org/10.1063/1.5019944
"""
Base.@kwdef struct AR1Unidir{V, A, B, C, NX, NY, R} <: DiscreteDefinition
    xi::V = [0.2, 0.3]
    a₁::A = 0.90693
    b₁::B = 0.40693
    c_xy::C = 0.5
    nx::NX = Normal(0, 0.40662)
    ny::NY = Normal(0, 0.40662)
    rng::R = Random.default_rng()
end

function system(definition::AR1Unidir)
    return DiscreteDynamicalSystem(eom_ar1_unidir, definition.xi, definition)
end

function eom_ar1_unidir(u, p::AR1Unidir, t)
    (; xi, a₁, b₁, c_xy, nx, ny, rng) = p
    x, y = u
    dx = a₁*x + rand(rng, nx)
    dy = b₁*y + c_xy*x + rand(rng, ny)
    return SVector{2}(dx, dy)
end
