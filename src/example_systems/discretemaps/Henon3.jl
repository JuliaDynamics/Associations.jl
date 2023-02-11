using LabelledArrays: @LArray
using StaticArrays: SVector, MVector
using DynamicalSystemsBase: DiscreteDynamicalSystem
using StateSpaceSets: Dataset

export Henon3

"""
    Henon3() <: DiscreteDefinition
    Henon3(; a = 0.1, b = 0.3, c = 0.1, xi = [0.1, 0.2, 0.3])

`Henon3` is a [`DiscreteDefinition`](@ref) definition for a lagged discrete dynamical
system consisting of three coupled 1D Henon maps (Papana et al., 2013)[^Papana2013].

## Equations of motion

```math
\\begin{aligned}
x_1(t+1) &= a - x_1(t)^2 + b x_1(t-2) \\
x_2(t+1) &= a - c x_1(t) x_2(t)- (1 - c) x_2(t)^2 + b x_2(t-1) \\
x_3(t+1) &= c x_2(t) x_3(t) - (1 - c) x_3(t)^2 + b x_3(t-1)
\\end{aligned}
```

Here ``c`` is the coupling constant. The system becomes completely synchronized
for ``c >= 0.7``. The initial condition `xi` is repeated over the first two time steps
before iteration starts.

[^Papana2013]:
    Papana, A., Kyrtsou, C., Kugiumtzis, D., & Diks, C. (2013). Simulation study of
    direct causality measures in multivariate time series. Entropy, 15(7), 2635–2661.
"""
struct Henon3{T, S, A, B, C} <: DiscreteDefinition
    m₁::MVector{2, T} # holds past states of x1
    m₂::MVector{2, T} # holds past states of x2
    m₃::MVector{2, T} # holds past states of x3
    xi::S
    a::A
    b::B
    c::C

    function Henon3(; a::A = 1.4, b::B = 0.3, c::C = 0.1,
            xi::S = [0.4, 0.5, 0.6]) where {A, B, C, S}
        T = eltype(1.0)
        m₁ = MVector{2, T}(repeat([xi[1]], 2))
        m₂ = MVector{2, T}(repeat([xi[2]], 2))
        m₃ = MVector{2, T}(repeat([xi[3]], 2))
        return new{T, S, A, B, C}(m₁, m₂, m₃, xi, a, b, c)
    end
end

function system(definition::Henon3)
    return DiscreteDynamicalSystem(eom_henon3, definition.xi, definition)
end

function eom_henon3(u, p::Henon3, t)
    # `u` is simply ignored here, because the state is stored in the memory vectors
    m₁, m₂, m₃ = p.m₁, p.m₂, p.m₃
    x₁₁, x₁₂ = m₁[1], m₁[2]
    x₂₁, x₂₂ = m₂[1], m₂[2]
    x₃₁, x₃₂ = m₃[1], m₃[2]

    a, b, c = p.a, p.b, p.c
    dx₁= a - x₁₁^2 + b*x₁₂
    dx₂= a - c*x₁₁*x₂₁ - (1 - c)*x₂₁^2 + b*x₂₂
    dx₃= a - c*x₂₁*x₃₁ - (1 - c)*x₃₁^2 + b*x₃₂

    new_state = SVector{3}(dx₁, dx₂, dx₃)
    update_state!(p, new_state) # Update memory vectors
    return new_state
end

function update_state!(p::Henon3, xnew::SVector{3})
    p.m₁[2] = p.m₁[1]
    p.m₁[1] = xnew[1]
    p.m₂[2] = p.m₂[1]
    p.m₂[1] = xnew[2]
    p.m₃[2] = p.m₃[1]
    p.m₃[1] = xnew[3]
end
