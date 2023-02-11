import Distributions: Uniform
import StaticArrays: SVector
import DynamicalSystemsBase: DiscreteDynamicalSystem
using Random
export Ikeda

"""
    Ikeda <: DiscreteDefinition
    Ikeda(; xi = [0.19, 0.21], c_xy = 1.0, c_yx = 1.0, a = 0.8, b = 12, c = 0.9,
        r₁ = 0.2, r₂ = 0.15, σ = 0.05, rng = Random.default_rng())

Initialise a discrete two-dimensional Ikeda map system, adapted from Cao et al.
(1997)[^Cao1997], by adding a noise term and allowing the influences from ``x \\to y``
(``c_{xy}``) and from ``y \\to x`` (``c_{yx}``) to be adjusted.

The difference equations are

```math
\\begin{aligned}
x(t+1) = 1 + \\mu(x \\cos{(\\theta)} - c_{yx} y \\sin{(\\theta)}) - min(\\dfrac{\\sigma \\xi_{t}^{(1)})}{(1-x)}, \\xi_{t}^{(2)} \\\\
y(t+1) = \\mu(y \\cos{(\\theta)} - c_{xy} x \\sin{(\\theta)}) - min(\\dfrac{\\sigma \\zeta_{t}^{(1)})}{(1-y)}, \\zeta_{t}^{(2)}
\\end{aligned}
```

## References

[^Cao1997]:
    Cao, Liangyue, Alistair Mees, and Kevin Judd. "Modeling and predicting
    non-stationary time series." International Journal of Bifurcation and
    Chaos 7.08 (1997): 1823-1831.
"""
Base.@kwdef struct Ikeda{V, C1, C2, A, B, C, R1, R2, Σ, R} <: DiscreteDefinition
    xi::V = [0.19, 0.21]
    c_xy::C1 = 1.0
    c_yx::C2 = 1.0
    a::A = 0.8
    b::B = 12
    c::C = 0.9
    r₁::R1 = 0.2
    r₂::R2 = 0.15
    σ::Σ = 0.05
    rng::R = Random.default_rng()
end

function system(definition::Ikeda)
    return DiscreteDynamicalSystem(eom_ikeda, definition.xi, definition)
end

function eom_ikeda(u, p::Ikeda, t)
    x, y = u
    (; xi, c_xy, c_yx, a, b, c, r₁, r₂, σ, rng) = p
    θ = a - b/(c + x^2 + y^2)
    μ = r₁*sin(t) - r₂
    d = Uniform(0.1, 0.4)

    dx = 1 + μ*(x*cos(θ) - c_yx*y*sin(θ)) - min(σ*rand(rng, d)/(1-x), rand(rng, d))
    dy = μ*(y*cos(θ) + c_xy*x*sin(θ)) -  min(σ*rand(rng, d)/(1-y), rand(rng, d))

    return SVector{2}(dx, dy)
end
