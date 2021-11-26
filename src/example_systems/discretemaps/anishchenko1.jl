using LabelledArrays

export anishchenko1

"""
    eom_anishchenko1(u, p, t)

Equations of motion for the system defined by eq. 13 in [1], which can 
give strange, nonchaotic attractors.

## Equations of motion 

The equations of motion are 

```math
\\begin{aligned}
dx &= \\alpha (1-s \\cos (2 \\pi \\phi )) \\cdot x(1-x) \\\\
dϕ &= (\\phi + \\omega ) \\mod{1}
\\end{aligned}
```

## References

1. Anishchenko, Vadim S., and Galina I. Strelkova. "Irregular attractors."
    Discrete dynamics in Nature and Society 2.1 (1998): 53-72.
"""
function eom_anishchenko1(u, p, t)
    x, ϕ = (u...,)
    α, s, ω = (p...,)
    dx = α*(1 - s*cos(2*pi*ϕ))*x*(1 - x)
    dϕ = (ϕ + ω) % 1

    return SVector{2}(dx, dϕ)
end

function anishchenko1(u₀, α, s, ω)
    p = @LArray [α, s, ω] (:α, :s, :ω)
    DiscreteDynamicalSystem(eom_anishchenko1, u₀, p)
end

"""
    anishchenko1(;u₀ = rand(2), α =3.277, s=0.1, ω=0.5*(sqrt(5)-1)) → DiscreteDynamicalSystem

Initialise the system defined by eq. 13 in [1], which can give strange, 
nonchaotic attractors.

## Equations of motion 

The equations of motion are 

```math
\\begin{aligned}
dx &= \\alpha (1-s \\cos (2 \\pi \\phi )) \\cdot x(1-x) \\\\
dϕ &= (\\phi + \\omega ) \\mod{1}
\\end{aligned}
```

## References

1. Anishchenko, Vadim S., and Galina I. Strelkova. "Irregular attractors."
    Discrete dynamics in Nature and Society 2.1 (1998): 53-72.
"""
anishchenko1(;u₀ = rand(2), α =3.277, s=0.1, ω=0.5*(sqrt(5)-1)) =
    anishchenko1(u₀, α, s, ω)
