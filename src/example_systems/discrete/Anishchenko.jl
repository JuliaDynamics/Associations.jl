using LabelledArrays: @LArray
using StaticArrays: SVector
using DynamicalSystemsBase: DiscreteDynamicalSystem

export Anishchenko

"""
    Anishchenko <: DiscreteDefinition
    Anishchenko(;u₀ = rand(2), α =3.277, s=0.1, ω=0.5*(sqrt(5)-1)) → DiscreteDynamicalSystem

Initialise the system defined by eq. 13 in [Anishchenko1998](@cite),
which can give strange, nonchaotic attractors.

## Equations of motion

```math
\\begin{aligned}
dx &= \\alpha (1-s \\cos (2 \\pi \\phi )) \\cdot x(1-x) \\\\
dϕ &= (\\phi + \\omega ) \\mod{1}
\\end{aligned}
```
"""
Base.@kwdef struct Anishchenko{V, A, S, Ω} <: DiscreteDefinition
    xi::V = [0.2, 0.3]
    α::A = 3.277
    s::S = 0.1
    ω::Ω = 0.5 * (sqrt(5) - 1)
end

function system(definition::Anishchenko)
    return DiscreteDynamicalSystem(eom_anischenko, definition.xi, definition)
end

function eom_anischenko(u, p::Anishchenko, t)
    x, ϕ = u
    α, s, ω = p.α, p.s, p.ω
    dx = α * (1 - s * cos(2*pi*ϕ)) * x * (1 - x)
    dϕ = (ϕ + ω) % 1

    return SVector{2}(dx, dϕ)
end
