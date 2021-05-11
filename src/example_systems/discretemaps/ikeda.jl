import Distributions: Uniform
import StaticArrays: SVector
import DynamicalSystems: DiscreteDynamicalSystem
using LabelledArrays

export ikeda

"""
    eom_ikeda(u, p, t)

Equations of motion for a discrete two-dimensional Ikeda map system, adapted from [1]
by adding a noise term and allowing the influences from ``x \\to y`` (``c_{xy}``) and 
from ``y \\to x`` (``c_{yx}``) to be adjusted. The difference equations are

```math
\\begin{aligned}
x(t+1) = 1 + \\mu(x(t) \\cos{(\\theta)} - c_{yx} y(t) \\sin{(\\theta)}) - min(\\dfrac{\\sigma \\xi_{t}^{(1)})}{(1-x)}, \\xi_{t}^{(2)} \\\\
y(t+1) = \\mu(y(t) \\cos{(\\theta)} - c_{xy} x(t) \\sin{(\\theta)}) - min(\\dfrac{\\sigma \\zeta_{t}^{(1)})}{(1-y)}, \\zeta_{t}^{(2)}
\\end{aligned}
```

## References

1. Cao, Liangyue, Alistair Mees, and Kevin Judd. "Modeling and predicting 
    non-stationary time series." International Journal of Bifurcation and 
    Chaos 7.08 (1997): 1823-1831.
"""
function eom_ikeda(u, p, t)
    x, y = u[1], u[2]
    c_xy, c_yx, a, b, c, r₁, r₂, σ = p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8]
    
    θ = a - b/(c + x^2 + y^2)
    μ = r₁*sin(t) - r₂
    d = Uniform(0.1, 0.4)
    
    dx = 1 + μ*(x*cos(θ) - c_yx*y*sin(θ)) - min(σ*rand(d)/(1-x), rand(d))
    dy = μ*(y*cos(θ) + c_xy*x*sin(θ)) -  min(σ*rand(d)/(1-y), rand(d))
    
    SVector{2}(dx, dy)
end

"""
    ikeda(; u₀ = rand(2), c_xy = 1.0, c_yx = 1.0, a = 0.8, b = 12, c = 0.9,
        r₁ = rand(Uniform(0.01, 0.3)), r₂ = rand(Uniform(0.01, 0.3)), σ = 0.05)

Initialise a discrete two-dimensional Ikeda map system, adapted from [1]
by adding a noise term and allowing the influences from ``x \\to y`` (``c_{xy}``) and 
from ``y \\to x`` (``c_{yx}``) to be adjusted.

As a rule-of-thumb, if parameters `a`, `b`, and `c` are drawn from uniform 
distributions on `[0.8, 1.5]`, `[10, 14]` and `[0.1, 0.9]`.

The difference equations are

```math
\\begin{aligned}
x(t+1) = 1 + \\mu(x \\cos{(\\theta)} - c_{yx} y \\sin{(\\theta)}) - min(\\dfrac{\\sigma \\xi_{t}^{(1)})}{(1-x)}, \\xi_{t}^{(2)} \\\\
y(t+1) = \\mu(y \\cos{(\\theta)} - c_{xy} x \\sin{(\\theta)}) - min(\\dfrac{\\sigma \\zeta_{t}^{(1)})}{(1-y)}, \\zeta_{t}^{(2)}
\\end{aligned}
```

## References

1. Cao, Liangyue, Alistair Mees, and Kevin Judd. "Modeling and predicting 
    non-stationary time series." International Journal of Bifurcation and 
    Chaos 7.08 (1997): 1823-1831.
"""
function ikeda(u₀, c_xy, c_yx, a, b, c, r₁, r₂, σ)
    p = @LArray [c_xy, c_yx, a, b, c, r₁, r₂, σ] (:c_xy, :c_yx, :a, :b, :c, :r₁, :r₂, :σ)
    DiscreteDynamicalSystem(eom_ikeda, u₀, p)
end

function ikeda(; u₀ = rand(2), c_xy = 1.0, c_yx = 1.0, a = 0.8, b = 12, c = 0.9,
        r₁ = rand(Uniform(0.01, 0.3)), r₂ = rand(Uniform(0.01, 0.3)), σ = 0.05)
    ikeda(u₀, c_xy, c_yx, a, b, c, r₁, r₂, σ)
end

export eom_ikeda, ikeda