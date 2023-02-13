using LabelledArrays: @LArray
using StaticArrays: SVector
using DynamicalSystemsBase: DiscreteDynamicalSystem
using Random

export Logistic4Chain

"""
    Logistic4Chain <: DiscreteDefinition
    Logistic4Chain(; xi = rand(4),
        rx = 3.9, ry = 3.6, rz = 3.6, rw = 3.8,
        cxy = 0.4, cyz = 0.4, cyw = 0.35,
        rng = Random.default_rng())

A variant of [`Logistic2Bidir`](@ref) where four variables `X`, `Y`, `Z`, `W`
are coupled in a chain `X → Y → Z → W` with dynamical noise.

## Description

The equations of motion are

```math
\\begin{align*}
x(t+1) &= r_x x(t)(1 - x(t))  \\\\
y(t+1) &= r_y f_{xy}^{t}(1 - f_{xy}^{t}) \\\\
z(t+1) &= r_z f_{yz}^{t}(1 - f_{yz}^{t}) \\\\
w(t+1) &= r_w f_{zw}^{t}(1 - f_{zw}^{t}) \\\\
f_{xy}^t &= \\dfrac{y(t) + c_{xy}(x(t) + \\sigma_{xy} \\xi_{xy}^t )}{1 + c_{xy} (1 + \\sigma_{xy} )} \\\\
f_{yz}^t &= \\dfrac{z(t) + c_{yz}(y(t) + \\sigma_{yz} \\xi_{yz}^t )}{1 + c_{yz} (1 + \\sigma_{yz} )}, \\\\
f_{zw}^t &= \\dfrac{w(t) + c_{zw}(z(t) + \\sigma_{zw} \\xi_{zw}^t )}{1 + c_{zw} (1 + \\sigma_{zw} )},
\\end{align*}
```

where `c_{xy}`, `c_{yz}`, `c_{zw}` controls the coupling strength from `x` to `y`, `y` to
`z`, and `z` to `w`, respectively. The magnitude of dynamical noise in these couplings
are controlled by ``\\sigma_{xy}``, ``\\sigma_{yz}``, and ``\\sigma_{zw}``, respectively.
``\\xi_{xy}``, ``\\xi_{yz}``, and ``\\xi_{zw}`` are noise terms that each iteration
are drawn from independent uniform distributions over the unit interval.

With default parameters and all dynamical noise terms set to zero, this is the system
from Ye et al. (2015)[^Ye2015] (but note that for some initial conditions,
this system wanders off to ``\\pm \\infty`` for some of the variables).

[^Ye2015]:
    Ye, Hao, et al. "Distinguishing time-delayed causal interactions using
    convergent cross mapping." Scientific reports 5 (2015): 14750
"""
Base.@kwdef struct Logistic4Chain{V, RX, RY, RZ, RW, C1, C2, C3, Σ1, Σ2, Σ3, RNG} <: DiscreteDefinition
    xi::V = [0.1, 0.2, 0.3, 0.4]
    rx::RX = 3.9
    ry::RY = 3.6
    rz::RZ = 3.6
    rw::RW = 3.8
    c_xy::C1 = 0.4
    c_yz::C2 = 0.4
    c_zw::C3 = 0.35
    σ_xy::Σ1 = 0.05
    σ_yz::Σ2 = 0.05
    σ_zw::Σ3 = 0.05
    rng::RNG = Random.default_rng()
end

function system(definition::Logistic4Chain)
    return DiscreteDynamicalSystem(eom_logistic4_chain, definition.xi, definition)
end

function eom_logistic4_chain(u, p::Logistic4Chain, t)
    (; xi, rx, ry, rz, rw, c_xy, c_yz, c_zw, σ_xy, σ_yz, σ_zw, rng) = p
    x, y, z, w = u
    f_xy = (y +  c_xy*(x + σ_xy * rand(rng)) ) / (1 + c_xy*(1+σ_xy))
    f_yz = (z +  c_yz*(y + σ_yz * rand(rng)) ) / (1 + c_yz*(1+σ_yz))
    f_zw = (w +  c_zw*(z + σ_zw * rand(rng)) ) / (1 + c_zw*(1+σ_zw))
    dx = rx * x * (1 - x)
    dy = ry * (f_xy) * (1 - f_xy)
    dz = rz * (f_yz) * (1 - f_yz)
    dw = rw * (f_zw) * (1 - f_zw)
    return SVector{4}(dx, dy, dz, dw)
end
