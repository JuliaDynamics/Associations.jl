using StaticArrays: SVector
using DynamicalSystemsBase: DiscreteDynamicalSystem
export Henon2

"""
    Henon2() <: DiscreteDefinition
    Henon2(;u₀ = [0.1, 0.2, 0.2, 0.3], c_xy = 2.0)

A bivariate system consisting of two identical 1D Henon maps with
unidirectional forcing ``X \\to Y `` (Krakovská et al., 2018)[^Krakovská2018].

## Equations of motion

The equations of motion are

```math
\\begin{align*}
x_1(t+1) &= 1.4 - x_1^2(t) + 0.3x_2(t) \\\\
x_2(t+1) &= x_1(t) \\\\
y_1(t+1) &= 1.4 - [c_{xy} x_1(t) y_1(t) + (1-c_{xy}) y_1^2(t)] + 0.3 y_2(t) \\\\
y_2(t+1) &= y_1(t)
\\end{align*}
```

[^Krakovská2018]:
    Krakovská, A., Jakubík, J., Chvosteková, M., Coufal, D., Jajcay, N., & Paluš, M. (2018).
    Comparison of six methods for the detection of causality in a bivariate time series.
    Physical Review E, 97(4), 042207.
"""
Base.@kwdef struct Henon2{R, C, V} <: DiscreteDefinition
    r::R = 3.4
    c_xy::C = 1.4
    xi::V = [0.1, 0.2, 0.2, 0.3]
end

function system(definition::Henon2)
    return DiscreteDynamicalSystem(eom_henon2, definition.xi, definition)
end

function eom_henon2(x, p::Henon2, t)
    c_xy = p.c_xy
    x₁, x₂, y₁, y₂ = x
    dx₁ = 1.4 - x₁^2 + 0.3*x₂
    dx₂ = x₁
    dy₁ = 1.4 - (c_xy * x₁ * y₁  +  (1 - c_xy)*y₁^2) + 0.3*y₂
    dy₂ = y₁
    return SVector{4}(dx₁, dx₂, dy₁, dy₂)
end
