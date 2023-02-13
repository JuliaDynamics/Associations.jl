using StaticArrays: SVector
using DynamicalSystemsBase: DiscreteDynamicalSystem
using Distributions: Normal

export AR1Bidir

"""
    AR1Bidir <: DiscreteDefinition
    AR1Bidir(;xi = [0.2, 0.3], a₁ = 0.5, b₁ = 0.7, c_xy = 0.1, c_yx = 0.2,
        nx = Normal(0, 0.3), ny = Normal(0, 0.3),
        rng::R = Random.default_rng())

A system consisting of two mutually coupled first order autoregressive processes.

## Equations of motion

```math
\\begin{aligned}
x(t+1) &= a_{1}x + c_{yx}y + \\epsilon_{x} \\\\
y(t+1) &= b_{1}y + c_{xy}x + \\epsilon_{y}
\\end{aligned}
```

where at each time step, ``\\epsilon_{x}`` and ``\\epsilon_{y}`` are drawn
from independent normal distributions `nx` and `ny`, respectively.
"""
Base.@kwdef struct AR1Bidir{V, A, B, C1, C2, NX, NY, R} <: DiscreteDefinition
    xi::V = [0.2, 0.3]
    a₁::A = 0.5
    b₁::B = 0.7
    c_xy::C1 = 0.1
    c_yx::C2 = 0.2
    nx::NX = Normal(0, 0.3)
    ny::NY = Normal(0, 0.3)
    rng::R = Random.default_rng()
end

function system(definition::AR1Bidir)
    return DiscreteDynamicalSystem(eom_ar1_bidir, definition.xi, definition)
end

function eom_ar1_bidir(u, p::AR1Bidir, t)
    (; xi, a₁, b₁, c_xy, c_yx, nx, ny, rng) = p
    x, y = u
    dx = a₁*x + c_yx*y + rand(rng, nx)
    dy = b₁*y + c_xy*x + rand(rng, ny)
    return SVector{2}(dx, dy)
end
