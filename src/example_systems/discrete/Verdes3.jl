using DynamicalSystemsBase: DiscreteDynamicalSystem
using Distributions: Normal

export Verdes3

"""
    Verdes3 <: DiscreteDefinition
    Verdes3(; ui = [0.1, 0.15, 0.2], ωy = 315, ωz = 80, σx = 0.0, σy = 0.0, σz = 0.0)

A 3D system where the response X is a highly nonlinear combination
of Y and Z (Verdes, 2005)[^Verdes2005]. The forcings Y and Z involve sines and cosines, and
have different periods, which controlled by `ωy` and `ωz`.

```math
\\begin{align*}
x(t+1) &= \\dfrac{y(t)(18y(t) - 27y(t)^2 + 10)}{2} + z(t)(1-z(t)) + ηx \\\\
y(t+1) &= \\dfrac{(1 - \\dfrac{\\cos(2\\pi)}{\\omega y}t)}{2} + ηy \\\\
z(t+1) &= \\dfrac{(1 - \\dfrac{\\sin(2\\pi)}{\\omega z}t)}{2} + ηz
\\end{align*}
```
where ηx, ηy, ηz is gaussian noise with mean 0 and standard deviation `σx`, `σy`
and `σz`.

[^Verdes2005]:
    Verdes, P. F. "Assessing causality from multivariate time series." Physical
    Review E 72.2 (2005): 026222.
"""
Base.@kwdef struct Verdes3{V, Ωy, Ωz, N1, N2, N3, R} <: DiscreteDefinition
    xi::V = [0.1, 0.15, 0.2]
    ωy::Ωy = 315
    ωz::Ωz = 80
    σx::N1 = Normal(0, 0.01)
    σy::N2 = Normal(0, 0.01)
    σz::N3 = Normal(0, 0.01)
    rng::R = Random.default_rng()
end

function system(definition::Verdes3)
    return DiscreteDynamicalSystem(eom_verdes3, definition.xi, definition)
end

function eom_verdes3(u, p::Verdes3, t)
    x, y, z = u
    (; xi, ωy, ωz, σx, σy, σz, rng) = p

    dx = y*(18y - 27y^2 + 10)/2 + z*(1-z) + rand(rng, σx)
    dy = (1 - cos((2*pi/ωy) * t))/2 + rand(rng, σy)
    dz = (1 - sin((2*pi/ωz) * t))/2 + rand(rng, σz)
    return SVector{3}(dx, dy, dz)
end
