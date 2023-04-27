using LabelledArrays: @LArray
using StaticArrays: SVector
using DynamicalSystemsBase: DiscreteDynamicalSystem


export Logistic3CommonDriver

"""
    Logistic3CommonDriver() <: DiscreteDefinition
    Logistic3CommonDriver(; u₀ = [0.1, 0.2, 0.3],
        r = 4.0, σx = 0.05, σy = 0.05, σz = 0.05,
        rng = Random.default_rng())

A discrete dynamical system consisting of three coupled 1D logistic maps
representing the response of two independent dynamical variables to the
forcing from a common driver (Runge, 2018)[^Runge2018].
The dynamical influence goes in the directions ``Z \\to X`` and ``Z \\to Y``.

## Equations of motion

The equations of motion are

```math
\\begin{align*}
x(t+1) &= (x(t)(r - r x(t) - z(t) + \\sigma_x \\eta_x)) \\mod 1 \\\\
y(t+1) &= (y(t)(r - r y(t) - z(t) + \\sigma_y \\eta_y)) \\mod 1 \\\\
z(t+1) &= (z(t)(r - r z(t) + \\sigma_z \\eta_z)) \\mod 1
\\end{align*}
```

Dynamical noise may be added to each of the dynamical variables by tuning the
parameters `σz`, `σx` and `σz`. Default values for the parameters
`r₁`, `r₂` and `r₃` are set such that the system exhibits chaotic behaviour,
with `r₁ = r₂ = r₃ = 4`.

[^Runge2018]:
    Runge, Jakob. Causal network reconstruction from time series: From theoretical
    assumptions to practical estimation, Chaos 28, 075310 (2018);
    doi: 10.1063/1.5025050
"""
Base.@kwdef struct Logistic3CommonDriver{V, R, Σx, Σy, Σz, RNG} <: DiscreteDefinition
    xi::V = [0.1, 0.2, 0.3]
    r::R = 4.0
    σx::Σx = 0.05
    σy::Σy = 0.05
    σz::Σz = 0.05
    rng::RNG = Random.default_rng()
end

function system(definition::Logistic3CommonDriver)
    return DiscreteDynamicalSystem(eom_logistic3_commondriver, definition.xi, definition)
end

# Note: Until the `eom_logistic2_bidir` function is deprecated, this function must
# be called something different; otherwise the DiscreteDynamicalSystem constructor
# doesn't work.
function eom_logistic3_commondriver(u, p::Logistic3CommonDriver, t)
    (; xi, r, σx, σy, σz, rng) = p
    x, y, z = u
    ηx = rand(rng)
    ηy = rand(rng)
    ηz = rand(rng)
    dx = (x*(r - r*x - z + σx*ηx)) % 1
    dy = (y*(r - r*y - z + σy*ηy)) % 1
    dz = (z*(r - r*z + σz*ηz)) % 1
    return SVector{3}(dx, dy, dz)
end
