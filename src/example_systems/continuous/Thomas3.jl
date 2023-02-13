using DynamicalSystemsBase: ContinuousDynamicalSystem
using StaticArrays: SVector

export Thomas3

"""
    Thomas3 <: ContinuousDefinition
    Thomas3(; xi = [0.11, 0.09, 0.10], b = 0.20)

[Thomas' cyclically symmetric attractor](https://en.wikipedia.org/wiki/Thomas%27_cyclically_symmetric_attractor)
is a continuous dynamical system with three variables. It has a single free parameter
`b`, for which interesting behaviour occurs when `b ∈ (0, 1)`. In particular,
the system is chaotic whenever `b < 0.20`.

## Definition

```math
\\begin{align*}
\\dfrac{dx}{dt} = sin(y) - bx \\\\
\\dfrac{dy}{dt} = sin(z) - by \\\\
\\dfrac{dz}{dt} = sin(x) - bz
\\end{align}
```
"""
Base.@kwdef struct Thomas3{V, B} <: ContinuousDefinition
    xi::V = [0.11, 0.09, 0.10]
    b::B = 0.20
end

function system(definition::Thomas3)
    return ContinuousDynamicalSystem(eom_thomas3, definition.xi, definition)
end

function eom_thomas3(u, p::Thomas3, t)
    b = p.b
    x, y, z = u
    dx = sin(y) - b*x
    dy = sin(z) - b*y
    dz = sin(x) - b*z
    return SVector{3}(dx, dy, dz)
end
