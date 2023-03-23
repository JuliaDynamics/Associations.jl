using LabelledArrays: @LArray
using StaticArrays: SVector
using DynamicalSystemsBase: DiscreteDynamicalSystem
using Random

export Logistic2Bidir

"""
    Logistic2Bidir() <: DiscreteDefinition
    Logistic2Bidir(; ui = [0.5, 0.5], c_xy = 0.1, c_yx = 0.1, r₁ = 3.78, r₂ = 3.66,
        σ_xy = 0.05, σ_yx = 0.05,
        rng = Random.default_rng())

A bivariate system consisting of two 1D noisy logistic maps which are bidirectionally
interacting (Diego et al., 2019)[^Diego2019].

## Equations of motion

```math
\\begin{align*}
x(t+1) &= r_1 f_{yx}^{t}(1 - f_{yx}^{t}) \\\\
y(t+1) &= r_2 f_{xy}^{t}(1 - f_{xy}^{t}) \\\\
f_{xy}^t &= \\dfrac{y(t) + c_{xy}(x(t) + \\sigma_{xy} \\xi_{xy}^t )}{1 + c_{xy} (1 + \\sigma_{xy} )} \\\\
f_{yx}^t &= \\dfrac{x(t) + c_{yx}(y(t) + \\sigma_{yx} \\xi_{yx}^t )}{1 + c_{yx} (1 + \\sigma_{yx} )},
\\end{align*}
```

Here, the coupling strength ``c_{xy}`` controls how strongly species ``x`` influences species
``y``, and vice versa for ``c_{yx}``. To simulate time-varying influence of unobserved
processes, we use the dynamical noise terms ``\\xi_{xy}^t`` and ``\\xi_{yx}^t``, drawn from a
uniform distribution with support on ``[0, 1]``. If ``\\sigma_{xy} > 0``, then the influence
of ``x`` on ``y`` is masked by dynamical noise equivalent to ``\\sigma_{xy} \\xi_{xy}^{t}`` at
the ``t``-th iteration of the map, and vice versa for ``\\sigma_{yx}``.

[^Diego2019]:
    Diego, David, Kristian Agasøster Haaga, and Bjarte Hannisdal. "Transfer entropy
    computation using the Perron-Frobenius operator." Physical Review E 99.4 (2019): 042212.
"""
Base.@kwdef struct Logistic2Bidir{V, C1, C2, R1, R2, Σx, Σy, R} <: DiscreteDefinition
    xi::V = [0.5, 0.5]
    c_xy::C1 = 0.1
    c_yx::C2 = 0.1
    r₁::R1 = 3.78
    r₂::R2 = 3.66
    σ_xy::Σx = 0.05
    σ_yx::Σy = 0.05
    rng::R = Random.default_rng()
end

function system(definition::Logistic2Bidir)
    return DiscreteDynamicalSystem(eom_logistic2bidir, definition.xi, definition)
end

# Note: Until the `eom_logistic2_bidir` function is deprecated, this function must
# be called something different; otherwise the DiscreteDynamicalSystem constructor
# doesn't work.
function eom_logistic2bidir(u, p::Logistic2Bidir, t)
    (; xi, c_xy, c_yx, r₁, r₂, σ_xy, σ_yx, rng) = p
    x, y = u
    f_xy = (y +  c_xy*(x + σ_xy * rand(rng)) ) / (1 + c_xy*(1+σ_xy))
    f_yx = (x +  c_yx*(y + σ_yx * rand(rng)) ) / (1 + c_yx*(1+σ_yx))
    dx = r₁ * (f_yx) * (1 - f_yx)
    dy = r₂ * (f_xy) * (1 - f_xy)
    return SVector{2}(dx, dy)
end


Base.@kwdef struct Logistic3CommonCause{V, C1, C2, R1, R2, R3, Σx, Σy, R} <: DiscreteDefinition
    xi::V = [0.5, 0.5, 0.5]
    c_zx::C1 = 0.1
    c_zy::C2 = 0.1
    rx::R1 = 4.0
    ry::R2 = 4.0
    rz::R3 = 4.0
    σ_zx::Σy = 0.05
    σ_zy::Σx = 0.05
    rng::R = Random.default_rng()
end

function system(definition::Logistic3CommonCause)
    return DiscreteDynamicalSystem(eom_logistic3_commoncause, definition.xi, definition)
end

function eom_logistic3_commoncause(u, p::Logistic3CommonCause, t)
    (; xi, c_zx, c_zy, rx, ry, rz, σ_zx, σ_zy, rng) = p
    x, y, z = u
    f_zx = (x +  c_zx*(z + σ_zx * rand(rng)) ) / (1 + c_zx*(1+σ_zx))
    f_zy = (y +  c_zy*(z + σ_zy * rand(rng)) ) / (1 + c_zy*(1+σ_zy))
    dx = rx * (f_zx) * (1 - f_zx)
    dy = ry * (f_zy) * (1 - f_zy)
    dz = rz * z * (1 - z)
    return SVector{3}(dx, dy, dz)
end

Base.@kwdef struct Logistic3CommonForcing{V, C1, C2, R1, R2, R3, Σx, Σy, R} <: DiscreteDefinition
    xi::V = [0.5, 0.5, 0.5]
    c_xz::C1 = 0.3
    c_yz::C2 = 0.3
    rx::R1 = 4.0
    ry::R2 = 4.0
    rz::R3 = 4.0
    σ_xz::Σy = 0.05
    σ_yz::Σx = 0.05
    rng::R = Random.default_rng()
end

function system(definition::Logistic3CommonForcing)
    return DiscreteDynamicalSystem(eom_logistic3_commonforcing, definition.xi, definition)
end

function eom_logistic3_commonforcing(u, p::Logistic3CommonForcing, t)
    (; xi, c_xz, c_yz, rx, ry, rz, σ_xz, σ_yz, rng) = p
    x, y, z = u
    f_xy_z = (z +  c_xz*(x + σ_xz * rand(rng)) + c_yz*(y + σ_yz * rand(rng))) / (1 + c_xz*(1+σ_yz) + c_yz*(1+σ_yz))
    dx = rx * x * (1 - x)
    dy = ry * y * (1 - y)
    dz = rz * (f_xy_z) * (1 - f_xy_z)
    return SVector{3}(dx, dy, dz)
end

export Logistic3CommonForcing
export Logistic3BidirExogenous

Base.@kwdef struct Logistic3BidirExogenous{V, C1, C2, R1, R2, R3, Σx, Σy, R} <: DiscreteDefinition
    xi::V = [0.5, 0.5]
    c_xy::C1 = 0.1
    c_yx::C2 = 0.1
    rx::R1 = 4.0
    ry::R2 = 4.0
    rz::R3 = 4.0
    σ_xy::Σx = 0.05
    σ_yx::Σy = 0.05
    rng::R = Random.default_rng()
end

function system(definition::Logistic3BidirExogenous)
    return DiscreteDynamicalSystem(eom_logistic3bidirexo, definition.xi, definition)
end

function eom_logistic3bidirexo(u, p::Logistic3BidirExogenous, t)
    (; xi, c_xy, c_yx, rx, ry, rz, σ_xy, σ_yx, rng) = p
    x, y, z = u
    f_xy = (y +  c_xy*(x + σ_xy * rand(rng)) ) / (1 + c_xy*(1+σ_xy))
    f_yx = (x +  c_yx*(y + σ_yx * rand(rng)) ) / (1 + c_yx*(1+σ_yx))
    dx = rx * (f_yx) * (1 - f_yx)
    dy = ry * (f_xy) * (1 - f_xy)
    dz = rz
    return SVector{3}(dx, dy, dz)
end
