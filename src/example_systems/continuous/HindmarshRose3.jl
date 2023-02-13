using StaticArrays: SVector
using DynamicalSystemsBase: ContinuousDynamicalSystem

export HindmarshRose3

"""
    HindmarshRose3 <: ContinuousDefinition
    HindmarshRose3(; xi = [0.1, 0.2, 0.3, p)

Initialise a Hindmarsh-Rose system, which is a model of neuronal
spiking.

## Description

```math
\\begin{align*}
\\dfrac{dx}{dt} &= y + \\phi(x) - z + I \\\\
\\dfrac{dy}{dt} &= \\psi(x) - y \\\\
\\dfrac{dz}{dt} &= r[s(x - x_R) - z],
\\end{align*}
```
where

```math
\\begin{aligned}
\\phi(x) &= -ax^3+bx^2
\\psi(x) &= c - dx^2
\\end{aligned}
```

If parameters other than the defaults are to be used, they must be
provided as a vector `[a, b, c, d, r, s, xᵣ, I]`.
"""
Base.@kwdef struct HindmarshRose3{V,A,B,C,D,R,S,X,II} <: ContinuousDefinition
    xi::V = [0.1, 0.2, 0.3]
    a::A = 1
    b::B = 3
    c::C = 1
    d::D = 5
    r::R = 1e-3
    s::S = 4
    xᵣ::X = -8/5
    I::II = -8
end

function system(definition::HindmarshRose3)
    return ContinuousDynamicalSystem(eom_hindmarshrose, definition.xi, definition)
end

function eom_hindmarshrose(u, p::HindmarshRose3, t)
    (; xi, a, b, c, d, r, s, xᵣ, I) = p
    x, y, z = u

	ϕ = -a*x^3 + b*x^2
	ψ = c - d*x^2
    dx = y + ϕ - z + I
	dy = ψ - y
	dz = r*(s*(x - xᵣ) - z)
    return SVector{3}(dx, dy, dz)
end
