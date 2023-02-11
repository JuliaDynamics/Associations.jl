using StaticArrays: SVector
using DynamicalSystemsBase: DiscreteDynamicalSystem
export ChaoticMaps3
"""
    ChaoticMaps3() <: DiscreteSystem
    ChaoticMaps3(; ui = [0.2, 0.1, 0.3], r = 3.4, c_xy = 0.5, c_xz = 0.5, c_yz = 0.3)

A model consisting of three coupled 1D maps, where ``x \\to y`` and ``x \\to z`` (Chen et
al., 2004)[^Chen2004].

## Equations of motion

```math
\\begin{aligned}
x(t) &= r x(t-1)( 1 - x(t-1)^2 ) e^{-x(t-1)^2} \\\\
y(t) &= r y(t-1)( 1 - y(t-1)^2 ) e^{-y(t-1)^2} + c_{xy} x(t-1) \\\\
z(t) &= r z(t-1)( 1 - z(t-1)^2 ) e^{-z(t-1)^2} + c_{xz} x(t-1) + c_{yz} y(t-1)
\\end{aligned}
```

The parameters `r`, `c_xy` and `c_yz` do not appear in the original paper,
but are added here for explorative purposes.

[^Chen2004]:
    Chen, Yonghong, et al. "Analyzing multiple nonlinear time series with extended Granger
    causality." Physics Letters A 324.1 (2004): 26-35
"""
Base.@kwdef struct ChaoticMaps3{R, V, C1 ,C2, C3} <: DiscreteSystem
    r::R = 3.4
    c_xy::C1 = 1.4
    c_xz::C2 = 0.3
    c_yz::C3 = 0.1
    xi::V = [0.2, 0.1, 0.3]
end

function system(definition::ChaoticMaps3)
    return DiscreteDynamicalSystem(eom_chaoticmaps3, definition.xi, definition)
end

function eom_chaoticmaps3(x, p::ChaoticMaps3, t)
    r, c_xy, c_xz, c_yz = p.r, p.c_xy, p.c_xz, p.c_yz
    x, y, z = x
    dx = r * x * (1 - x^2) * exp(-x^2)
    dy = r * y * (1 - y^2) * exp(-y^2) + c_xy * x
    dz = r * z * (1 - z^2) * exp(-z^2) + c_xz * x + c_yz * y
    return SVector{3}(dx, dy, dz)
end
