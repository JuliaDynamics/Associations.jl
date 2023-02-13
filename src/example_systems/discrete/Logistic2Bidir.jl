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
\\begin{align}
x(t+1) &= r_1 f_{yx}^{t}(1 - f_{yx}^{t}) \\\\
y(t+1) &= r_2 f_{xy}^{t}(1 - f_{xy}^{t}) \\\\
f_{xy}^t &= \\dfrac{y(t) + c_{xy}(x(t) + \\sigma_{xy} \\xi_{xy}^t )}{1 + c_{xy} (1 + \\sigma_{xy} )} \\\\
f_{yx}^t &= \\dfrac{x(t) + c_{yx}(y(t) + \\sigma_{yx} \\xi_{yx}^t )}{1 + c_{yx} (1 + \\sigma_{yx} )},
\\end{align}
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
