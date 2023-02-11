using LabelledArrays: @LArray
using StaticArrays: SVector
using DynamicalSystemsBase: DiscreteDynamicalSystem
using DynamicalSystemsBase: trajectory
using Distributions: Uniform

export Logistic2Unidir


"""
    Logistic2Unidir <: DiscreteDefinition
    Logistic2Unidir(;u₀ = rand(2), c_xy = 0.1, σ = 0.05,
        r₁ = 3.78, r₂ = 3.66)

A bivariate system consisting of two 1D noisy logistic maps which are undirectionally
coupled `x → y` (Diego et al., 2019)[^Diego2019].

## Equations of motion

The equations of motion are

```math
\\begin{aligned}
x(t+1) &= r_1 x(t)(1 - x(t)) \\\\
y(t+1) &= r_2 f(x,y)(1 - f(x,y)),
\\end{aligned}
```

with

```math
\\begin{aligned}
f(x,y) = \\dfrac{y + \\frac{c_{xy}(x \\xi )}{2}}{1 + \\frac{c_{xy}}{2}(1+ \\sigma )}
\\end{aligned}
```

The parameter `c_xy` controls how strong the dynamical forcing is. If `σ > 0`,
dynamical noise masking the influence of  `x` on `y` equivalent to
``\\sigma \\cdot \\xi`` is added at each iteration. Here,``\\xi`` is a draw from a
flat distribution on ``[0, 1]``. Thus, setting `σ = 0.05` is equivalent to
add dynamical noise corresponding to a maximum of ``5 \\%`` of the possible
range of values of the logistic map.

[^Diego2019]:
    Diego, David, Kristian Agasøster Haaga, and Bjarte Hannisdal. "Transfer entropy
    computation using the Perron-Frobenius operator." Physical Review E 99.4 (2019): 042212.
"""
Base.@kwdef struct Logistic2Unidir{V, C, R1, R2, Σy, R} <: DiscreteDefinition
    xi::V = [0.5, 0.5]
    c_xy::C = 0.1
    r₁::R1 = 3.78
    r₂::R2 = 3.66
    σ_xy::Σy = 0.05
    rng::R = Random.default_rng()
end

function system(definition::Logistic2Unidir)
    return DiscreteDynamicalSystem(eom_logistic2uni, definition.xi, definition)
end

# Note: Until the `eom_logistic2_bidir` function is deprecated, this function must
# be called something different; otherwise the DiscreteDynamicalSystem constructor
# doesn't work.
function eom_logistic2uni(u, p::Logistic2Unidir, t)
    (; xi, c_xy, r₁, r₂, σ_xy, rng) = p
    x, y = u
    f_xy = (y +  (c_xy*(x + σ_xy * rand(rng))/2) ) / (1 + (c_xy/2)*(1+σ_xy))

    dx = r₁ * x * (1 - x)
    dy = r₂ * (f_xy) * (1 - f_xy)
    return SVector{2}(dx, dy)
end
