using DynamicalSystemsBase: DiscreteDynamicalSystem

export UlamLattice

"""
    UlamLattice <: DiscreteDefinition
    UlamLattice(; D::Int = 10; ui = sin.(1:10), ε::Real = 0.10)

A lattice of `D` unidirectionally coupled ulam maps[^Schreiber2000] defined as

```math
x^{m}_{t+1} = f(\\epsilon x^{m-1}_{t} + (1 - \\epsilon) x_{t}^{m}),
```

where ``m = 1, 2, \\ldots, D`` and ``f(x) = 2 - x^2``. In this system, information transfer
happens only in the direction of increasing ``m``.

[^Schreiber2000]:
    Schreiber, Thomas. "Measuring information transfer." Physical review letters 85.2
    (2000): 461.
"""
Base.@kwdef struct UlamLattice{V, E, F} <: DiscreteDefinition
    xi::V = sin.(1:10)
    ε::E = 0.1
    f::F = x -> 2 - x^2
end

function system(definition::UlamLattice)
    return DiscreteDynamicalSystem(eom_ulamlattice, definition.xi, definition)
end

function eom_ulamlattice(dx, x, p::UlamLattice, t)
    (; xi, ε, f) = p
    # `u` is simply ignored here, because the state is stored in the memory vectors
    dx[1] = f(ε*x[length(dx)] + (1-ε)*x[1])
    for i in 2:length(dx)
        dx[i] = f(ε*x[i-1] + (1-ε)*x[i])
    end
end
