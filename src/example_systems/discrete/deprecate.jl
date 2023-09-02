
using LabelledArrays: @LArray
using StaticArrays: SVector
using DynamicalSystemsBase: DiscreteDynamicalSystem
using DynamicalSystemsBase: trajectory
using StateSpaceSets: StateSpaceSet
using Distributions: Normal, Uniform

export ar1_unidir
export henon_triple
export henon2
export ar1_bidir
export ikeda
export linearmap1
export logistic2_unidir
export logistic2_bidir
export logistic3
export logistic4
export nonlinear_3d
export nontrivial_pegiun
export ulam
export var1
export verdes

function eom_ar1_unidir(x, p, n)
    a₁, b₁, c_xy, σ = (p...,)
    x, y = (x...,)
    ξ₁ = rand(Normal(0, σ))
    ξ₂ = rand(Normal(0, σ))

    dx = a₁*x + ξ₁
    dy = b₁*y + c_xy*x + ξ₂
    return SVector{2}(dx, dy)
end

function ar1_unidir(u₀, a₁, b₁, c_xy, σ)
    @warn "`ar1_unidir` is deprecated in CausalityTools v2. "*
        "Use `system(AR1Unidir())` instead, which returns a `DiscreteDynamicalSystem` "*
        "that can be iterated"
    p = @LArray [a₁, b₁, c_xy, σ] (:a₁, :b₁, :c_xy, :σ)
    DiscreteDynamicalSystem(eom_ar1_unidir, u₀, p)
end

"""
    ar1_unidir(u₀, a₁ = 0.90693, b₁ = 0.40693, c_xy = 0.5,
        σ = 0.40662) → DiscreteDynamicalSystem

A bivariate, order one autoregressive model, where ``x \\to y`` (Paluš et al,
2018)[^Paluš2018].

## Equations of motion

```math
\\begin{aligned}
x(t+1) &= a_1 x(t) + \\xi_{1} \\\\
y(t+1) &= b_1 y(t) - c_{xy} x + \\xi_{2},
\\end{aligned}
```

where ``\\xi_{1}`` and ``\\xi_{2}`` are drawn from normal distributions
with zero mean and standard deviation `σ` at each iteration.

[^Paluš2018]:
    Paluš, M., Krakovská, A., Jakubík, J., & Chvosteková, M. (2018). Causality,
    dynamical systems and the arrow of time. Chaos: An Interdisciplinary Journal of
    Nonlinear Science, 28(7), 075307. http://doi.org/10.1063/1.5019944
"""
ar1_unidir(;u₀ = rand(2), a₁ = 0.90693, b₁ = 0.40693, c_xy = 0.5, σ = 0.40662) =
    ar1_unidir(u₀, a₁, b₁, c_xy, σ)

function eom_henon_triple(u, p, n)
    O = zeros(Float64, n + 3, 3)
    x₁, x₂, x₃ = (u...,)
    a, b, c = (p...,)

    # Propagate initial condition to the three first time steps.
    for i = 1:3
        O[i, 1] = x₁
        O[i, 2] = x₂
        O[i, 3] = x₃
    end
    for i = 4:n+3
        x₁1 = O[i-1, 1]
        x₁2 = O[i-2, 1]
        x₂1 = O[i-1, 2]
        x₂2 = O[i-2, 2]
        x₃1 = O[i-1, 3]
        x₃2 = O[i-2, 3]

        x₁new = a - x₁1^2 + b*x₁2
        x₂new = a - c*x₁1*x₂1 - (1 - c)*x₂1^2 + b*x₂2
        x₃new = a - c*x₂1*x₃1 - (1 - c)*x₃1^2 + b*x₃2

        O[i, 1] = x₁new
        O[i, 2] = x₂new
        O[i, 3] = x₃new
    end

    return O[4:end, :]
end


function henon_triple(u₀, a, b, c, n::Int, n_transient::Int)
    @warn "`henon_triple` is deprecated in CausalityTools v2. "*
        "Use `system(HenonTriple())` instead, which returns a `DiscreteDynamicalSystem` "*
        "that can be iterated."
    p = @LArray [a, b, c] (:a, :b, :c)
    o = eom_henon_triple(u₀, p, n + n_transient)
    x, y, z = o[n_transient+1:end, 1], o[n_transient+1:end, 2], o[n_transient+1:end, 3]
    return StateSpaceSet(x, y, z)
end

"""
    henon_triple(x, p, n) → Function

Iterate a 3D discrete system consisting of coupled Henon maps where the coupling
is x1 → x2 → x3 [1]. This version allows for tweaking the parameters of the
equations.

The difference equations are:

```math
\\begin{aligned}
x_1(t+1) &= a - x_1(t)^2 + b x_1(t-2) \\
x_2(t+1) &= a - c x_1(t) x_2(t)- (1 - c) x_2(t)^2 + b x_2(t-1) \\
x_3(t+1) &= c x_2(t) x_3(t) - (1 - c) x_3(t)^2 + b x_3(t-1)
\\end{aligned}
```

Here ``c`` is the coupling constant. The system becomes completely synchronized
for ``c >= 0.7``.

# References
1. Papana, A., Kyrtsou, C., Kugiumtzis, D., & Diks, C. (2013). Simulation study of
direct causality measures in multivariate time series. Entropy, 15(7), 2635–2661.
"""
function henon_triple(;u₀ = rand(3), a = 1.4, b = 0.3, c = 0.0,  n::Int = 100, n_transient::Int = 100)
    henon_triple(u₀, a, b, c, n, n_transient)
end


export henon2

function eom_henon2(x, p, n)
    c_xy = p[1]
    x₁, x₂, y₁, y₂ = x
    dx₁ = 1.4 - x₁^2 + 0.3*x₂
    dx₂ = x₁
    dy₁ = 1.4 - (c_xy * x₁ * y₁  +  (1 - c_xy)*y₁^2) + 0.3*y₂
    dy₂ = y₁
    return SVector{4}(dx₁, dx₂, dy₁, dy₂)
end

function henon2(u₀, c_xy)
    @warn "`henon2` is deprecated in CausalityTools v2. "*
    "Use `system(Henon2())` instead, which returns a `DiscreteDynamicalSystem` "*
    "that can be iterated."
    p = (c_xy,)
    DiscreteDynamicalSystem(eom_henon2, u₀, p)
end

"""
    henon2(;u₀ = [0.1, 0.2, 0.2, 0.3], c_xy = 2.0) → DiscreteDynamicalSystem

A bivariate system consisting of two identical 1D Henon maps with
unidirectional forcing ``X \\to Y `` (Krakovská et al., 2018)[^Krakovská2018].

## Equations of motion

The equations of motion are

```math
\\begin{aligned}
x_1(t+1) &= 1.4 - x_1^2(t) + 0.3x_2(t) \\\\
x_2(t+1) &= x_1(t) \\\\
y_1(t+1) &= 1.4 - [c_{xy} x_1(t) y_1(t) + (1-c_{xy}) y_1^2(t)] + 0.3 y_2(t) \\\\
y_2(t+1) &= y_1(t)
\\end{aligned}
```

[^Krakovská2018]:
    Krakovská, A., Jakubík, J., Chvosteková, M., Coufal, D., Jajcay, N., & Paluš, M. (2018).
    Comparison of six methods for the detection of causality in a bivariate time series.
    Physical Review E, 97(4), 042207.
"""
henon2(;u₀ = [0.1, 0.2, 0.2, 0.3], c_xy = 2.0) = henon2(u₀, c_xy)


function eom_ar1_bidir(x, p, t)
    a₁, b₁, c_xy, c_yx, ϵx, ϵy = p
    x, y = (x...,)
    dx = a₁*x + c_yx*y + rand(ϵx)
    dy = b₁*y + c_xy*x + rand(ϵy)
    return SVector{2}(dx, dy)
end

function ar1_bidir(u₀,a₁, b₁, c_xy, c_yx, σx, σy)
    @warn "`ar1_bidir` is deprecated in CausalityTools v2. "*
    "Use `system(AR1Bidir())` instead, which returns a `DiscreteDynamicalSystem` "*
    "that can be iterated."
    ϵx = Normal(0, σx)
    ϵy = Normal(0, σy)
    p = (a₁, b₁, c_xy, c_yx, ϵx, ϵy)
    return DiscreteDynamicalSystem(eom_ar1_bidir, u₀, p)
end

"""
    ar1_bidir(;u₀ = rand(2), a₁ = 0.5, b₁ = 0.7, c_xy = 0.1, c_yx = 0.2,
        σx = 0.3, σy = 0.3) → DiscreteDynamicalSystem

A system consisting of two mutually coupled first order autoregressive processes.

## Equations of motion

```math
\\begin{aligned}
x(t+1) &= a_{1}x + c_{yx}y + \\epsilon_{x} \\\\
y(t+1) &= b_{1}y + c_{xy}x + \\epsilon_{y}
\\end{aligned}
```

where at each time step, ``\\epsilon_{x}`` and ``\\epsilon_{y}`` are drawn
from independent normal distributions with zero mean and standard deviations `σx` and `σy`,
respectively.
"""
ar1_bidir(;a₁ = 0.5, b₁ = 0.7, u₀ = rand(2), c_xy = 0.1, c_yx = 0.2, σx = 0.3, σy = 0.3) =
    ar1_bidir(u₀, a₁, b₁, c_xy, c_yx, σx, σy)

export anishchenko1

function eom_anishchenko1(u, p, t)
    x, ϕ = (u...,)
    α, s, ω = (p...,)
    dx = α*(1 - s*cos(2*pi*ϕ))*x*(1 - x)
    dϕ = (ϕ + ω) % 1

    return SVector{2}(dx, dϕ)
end

function anishchenko1(u₀, α, s, ω)
    @warn "`anishchenko1` is deprecated in CausalityTools v2. "*
    "Use `system(Anishchenko())` instead, which returns a `DiscreteDynamicalSystem` "*
    "that can be iterated."
    p = @LArray [α, s, ω] (:α, :s, :ω)
    DiscreteDynamicalSystem(eom_anishchenko1, u₀, p)
end

"""
    anishchenko1(;u₀ = rand(2), α =3.277, s=0.1, ω=0.5*(sqrt(5)-1)) → DiscreteDynamicalSystem

Initialise the system defined by eq. 13 in [Anishchenko1998](@cite),
which can give strange, nonchaotic attractors.

## Equations of motion

```math
\\begin{aligned}
dx &= \\alpha (1-s \\cos (2 \\pi \\phi )) \\cdot x(1-x) \\\\
dϕ &= (\\phi + \\omega ) \\mod{1}
\\end{aligned}
```
"""
anishchenko1(;u₀ = rand(2), α =3.277, s=0.1, ω=0.5*(sqrt(5)-1)) =
    anishchenko1(u₀, α, s, ω)


export ikeda

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
    @warn "`ikeda` is deprecated in CausalityTools v2. "*
    "Use `system(Ikeda())` instead, which returns a `DiscreteDynamicalSystem` "*
    "that can be iterated."
    p = @LArray [c_xy, c_yx, a, b, c, r₁, r₂, σ] (:c_xy, :c_yx, :a, :b, :c, :r₁, :r₂, :σ)
    DiscreteDynamicalSystem(eom_ikeda, u₀, p)
end

function ikeda(; u₀ = rand(2), c_xy = 1.0, c_yx = 1.0, a = 0.8, b = 12, c = 0.9,
        r₁ = rand(Uniform(0.01, 0.3)), r₂ = rand(Uniform(0.01, 0.3)), σ = 0.05)
    ikeda(u₀, c_xy, c_yx, a, b, c, r₁, r₂, σ)
end

function eom_linearmap1(x, p, t)
    c = p[1]
    x, y = (x...,)
    t = t + 3
    dx = 3.4*x*(t - 1)*(1 - x^2*(t - 1))*exp(-x^2*(t - 1)) + 0.8*x*(t - 2) + rand(Normal(0, 0.05))
    dy = 3.4*y*(t - 1)*(1 - y^2*(t - 1))*exp(-y^2*(t - 1)) + 0.5*y*(t - 2) + c*x*(t - 2) + rand(Normal(0, 0.05))
    return SVector{2}(dx, dy)
end

function linearmap1(u₀, c)
    @warn "`linearmap1` is deprecated in CausalityTools v2. "*
    "Use `system(ChaoticNoisyLinear2())` instead, which returns a `DiscreteDynamicalSystem` "*
    "that can be iterated."
    p = @LArray [c] (:c)
    DiscreteDynamicalSystem(eom_linearmap2, u₀, p)
end

"""
    linearmap1(;u₀ = [1, rand(2)], c = 0.5) → DiscreteDynamicalSystem

Linear map from [Chen2004](@citet).
"""
linearmap1(;u₀ = rand(2), c = 0.5) = linearmap1(u₀, c)


export logistic2_bidir

"""
    logistic2_bidir(u₀, c_xy, c_yx, r₁, r₂, σ_xy, σ_yx)

Equations of motion for a bidirectional logistic model for the chaotic
population dynamics of two interacting species. This system is from [1],
and is given by

```math
\\begin{align}
x(t+1) &= r_1 f_{yx}^{t}(1 - f_{yx}^{t}) \\
y(t+1) &= r_2 f_{xy}^{t}(1 - f_{xy}^{t}) \\
f_{xy}^t &= \\dfrac{y(t) + c_{xy}(x(t) + \\sigma_{xy} \\xi_{xy}^t )}{1 + c_{xy} (1 + \\sigma_{xy} )} \\
f_{yx}^t &= \\dfrac{x(t) + c_{yx}(y(t) + \\sigma_{yx} \\xi_{yx}^t )}{1 + c_{yx} (1 + \\sigma_{yx} )},
\\end{align}
```

where the coupling strength ``c_{xy}`` controls how strongly species ``x`` influences species
``y``, and vice versa for ``c_{yx}``. To simulate time-varying influence of unobserved
processes, we use the dynamical noise terms ``\\xi_{xy}^t`` and ``\\xi_{yx}^t``, drawn from a
uniform distribution with support on ``[0, 1]``. If ``\\sigma_{xy} > 0``, then the influence
of ``x`` on ``y`` is masked by dynamical noise equivalent to ``\\sigma_{xy} \\xi_{xy}^{t}`` at
the ``t``-th iteration of the map, and vice versa for ``\\sigma_{yx}``.

## References

1. Diego, David, Kristian Agasøster Haaga, and Bjarte Hannisdal. "Transfer entropy computation
    using the Perron-Frobenius operator." Physical Review E 99.4 (2019): 042212.
"""
function eom_logistic2_bidir(dx, x, p, n)

    # c_xy is the coupling from x to y
    # c_yx is the coupling from y to x
    # σ_yx is the dynamical noise from y to x
    # σ_xy is the dynamical noise from y to x
    c_xy, c_yx, r₁, r₂, σ_xy, σ_yx = (p...,)

    ξ₁ = rand() # random number from flat distribution on [0, 1]
    ξ₂ = rand() # random number from flat distribution on [0, 1]
    x, y = x[1], x[2]

    f_xy = (y +  c_xy*(x + σ_xy*ξ₁) ) / (1 + c_xy*(1+σ_xy))
    f_yx = (x +  c_yx*(y + σ_yx*ξ₂) ) / (1 + c_yx*(1+σ_yx))

    dx[1] = r₁ * (f_yx) * (1 - f_yx)
    dx[2] = r₂ * (f_xy) * (1 - f_xy)
    return
end

function logistic2_bidir(u₀, c_xy, c_yx, r₁, r₂, σ_xy, σ_yx)
    @warn "`logistic2_bidir` is deprecated in CausalityTools v2. "*
    "Use `system(Logistic2Bidir())` instead, which returns a "*
    "`DiscreteDynamicalSystem` that can be iterated."
    p = @LArray [c_xy, c_yx, r₁, r₂, σ_xy, σ_yx] (:c_xy, :c_yx, :r₁, :r₂, :σ_xy, :σ_yx)
    DiscreteDynamicalSystem(eom_logistic2_bidir, u₀, p)
end

"""
    logistic2_bidir(;u₀ = rand(2), c_xy = 0.1, c_yx = 0.1,
        r₁ = 3.78, r₂ = 3.66, σ_xy = 0.05, σ_yx = 0.05)

A bidirectional logistic model for the chaotic population dynamics of two interacting
species [1].

## Equations of motion

The equations of motion are

```math
\\begin{align}
x(t+1) &= r_1 f_{yx}^{t}(1 - f_{yx}^{t}) \\\\
y(t+1) &= r_2 f_{xy}^{t}(1 - f_{xy}^{t}) \\\\
f_{xy}^t &= \\dfrac{y(t) + c_{xy}(x(t) + \\sigma_{xy} \\xi_{xy}^t )}{1 + c_{xy} (1 + \\sigma_{xy} )} \\\\
f_{yx}^t &= \\dfrac{x(t) + c_{yx}(y(t) + \\sigma_{yx} \\xi_{yx}^t )}{1 + c_{yx} (1 + \\sigma_{yx} )},
\\end{align}
```

where the coupling strength ``c_{xy}`` controls how strongly species ``x`` influences species
``y``, and vice versa for ``c_{yx}``. To simulate time-varying influence of unobserved
processes, we use the dynamical noise terms ``\\xi_{xy}^t`` and ``\\xi_{yx}^t``, drawn from a
uniform distribution with support on ``[0, 1]``. If ``\\sigma_{xy} > 0``, then the influence
of ``x`` on ``y`` is masked by dynamical noise equivalent to ``\\sigma_{xy} \\xi_{xy}^{t}`` at
the ``t``-th iteration of the map, and vice versa for ``\\sigma_{yx}``.
"""
logistic2_bidir(;u₀ = rand(2), c_xy = 0.1, c_yx = 0.1,
    r₁ = 3.78, r₂ = 3.66, σ_xy = 0.05, σ_yx = 0.05) =
    logistic2_bidir(u₀, c_xy, c_yx, r₁, r₂, σ_xy, σ_yx)

    export logistic2_unidir

"""
    eom_logistic2(dx, x, p, n) → function

Equations of motions for a system consisting of two coupled logistic maps where
X unidirectionally influences Y [1].


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

1. Diego, David, Kristian Agasøster Haaga, and Bjarte Hannisdal. "Transfer entropy computation
    using the Perron-Frobenius operator." Physical Review E 99.4 (2019): 042212.
"""

function eom_logistic2_unidir(dx, x, p, n)
    c_xy, r₁, r₂, σ = (p...,)
    ξ = rand() # random number from flat distribution on [0, 1]
    x, y = x[1], x[2]
    f_xy = (y +  (c_xy*(x + σ*ξ)/2) ) / (1 + (c_xy/2)*(1+σ))

    dx[1] = r₁ * x * (1 - x)
    dx[2] = r₂ * (f_xy) * (1 - f_xy)
    return
end

function logistic2_unidir(u₀, c_xy, r₁, r₂, σ)
    @warn "`logistic2_unidir` is deprecated in CausalityTools v2. "*
    "Use `system(Logistic2Unidir())` instead, which returns a "*
    "`DiscreteDynamicalSystem` that can be iterated."
    p = @LArray [c_xy, r₁, r₂, σ] (:c_xy, :r₁, :r₂, :σ)
    DiscreteDynamicalSystem(eom_logistic2_unidir, u₀, p)
end

"""
    logistic2(;u₀ = rand(2), c_xy = 0.1, σ = 0.05,
        r₁ = 3.78, r₂ = 3.66) → DiscreteDynamicalSystem

Initialise a system consisting of two coupled logistic maps where X
unidirectionally influences Y. By default, the parameters `r₁` and `r₂` are set
to values yielding chaotic behaviour.

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

## References

1. Diego, David, Kristian Agasøster Haaga, and Bjarte Hannisdal. "Transfer entropy computation
    using the Perron-Frobenius operator." Physical Review E 99.4 (2019): 042212.
"""
logistic2_unidir(;u₀ = rand(2), c_xy = 0.1, r₁ = 3.78, r₂ = 3.66, σ = 0.05) =
    logistic2_unidir(u₀, c_xy, r₁, r₂, σ)

#To get chaotic realisation, check that the orbit doesn't settle to a few unique values
function good_logistic_unidir_trajectory(npts::Int;
        Ttr = 1000, dt = 1,
        c_xy = 0.5,
        Dr₁ = Uniform(3.6, 4.0),
        Dr₂ = Uniform(3.6, 4.0),
        σ = 0.0,
        n_maxtries = 300)

    n_tries = 0
    while n_tries <= n_maxtries
        s = logistic2_unidir(u₀ = rand(2),
            c_xy = c_xy,
            σ = σ,
            r₁ = rand(Dr₁),
            r₂ = rand(Dr₂))

        o = trajectory(s, npts * dt - 1, Ttr = Ttr, dt = dt)

        # Ensure there are not too many repeated values, so we don't have trivial behaviour

        if length(unique(o[:, 1])) > npts * 0.9 && length(unique(o[:, 2])) > npts * 0.9
            return o
        end

        n_tries += 1
    end
end


export logistic3

"""
    eom_logistic3(u, p, t)

Equations of motion for a system consisting of three coupled logistic map
representing the response of two independent dynamical variables to the
forcing from a common driver. The dynamical influence goes in the directions
Z → X and Z → Y.

## Equations of motion

The equations of motion are

```math
\\begin{aligned}
x(t+1) = (x(t)(r - r_1 x(t) - z(t) + σ_x η_x)) \\mod 1 \\\\
y(t+1) = (y(t)(r - r_2 y(t) - z(t) + σ_y η_y)) \\mod 1 \\\\
z(t+1) = (z(t)(r - r_3 z(t) + σ_z η_z)) \\mod 1
\\end{aligned}
```

Dynamical noise may be added to each of the dynamical variables by tuning the
parameters `σz`, `σx` and `σz`. Default values for the parameters
`r₁`, `r₂` and `r₃` are set such that the system exhibits chaotic behaviour,
with `r₁ = r₂ = r₃ = 4`.

## References

1. Runge, Jakob. Causal network reconstruction from time series: From theoretical
    assumptions to practical estimation, Chaos 28, 075310 (2018);
    doi: 10.1063/1.5025050
"""
function eom_logistic3(u, p, t)
    r₁, r₂, r₃, σx, σy, σz = (p...,)
    x, y, z = (u...,)

    # Independent dynamical noise for each variable.
    ηx = rand()
    ηy = rand()
    ηz = rand()

    dx = (x*(r₁ - r₁*x - z + σx*ηx)) % 1
    dy = (y*(r₂ - r₂*y - z + σy*ηy)) % 1
    dz = (z*(r₃ - r₃*z + σz*ηz)) % 1
    return SVector{3}(dx, dy, dz)
end

function logistic3(u₀, r₁, r₂, r₃, σx, σy, σz)
    @warn "`logistic3` is deprecated in CausalityTools v2. "*
    "Use `system(Logistic3CommonDriver())` instead, which returns a "*
    "`DiscreteDynamicalSystem` that can be iterated."
    p = @LArray [r₁, r₂, r₃, σx, σy, σz] (:r₁, :r₂, :r₃, :σx, :σy, :σz)
    DiscreteDynamicalSystem(eom_logistic3, u₀, p)
end

"""
    logistic3(;u₀ = rand(3), r = 4,
        σx = 0.05, σy = 0.05, σz = 0.05) → DiscreteDynamicalSystem

Initialise a dynamical system consisting of three coupled logistic map
representing the response of two independent dynamical variables to the
forcing from a common driver. The dynamical influence goes in the directions
``Z \\to X`` and ``Z \\to Y``.

## Equations of motion

The equations of motion are

```math
\\begin{aligned}
x(t+1) = (x(t)(r - r_1 x(t) - z(t) + σ_x η_x)) \\mod 1 \\\\
y(t+1) = (y(t)(r - r_2 y(t) - z(t) + σ_y η_y)) \\mod 1 \\\\
z(t+1) = (z(t)(r - r_3 z(t) + σ_z η_z)) \\mod 1
\\end{aligned}
```

Dynamical noise may be added to each of the dynamical variables by tuning the
parameters `σz`, `σx` and `σz`. Default values for the parameters
`r₁`, `r₂` and `r₃` are set such that the system exhibits chaotic behaviour,
with `r₁ = r₂ = r₃ = 4`.

## References

1. Runge, Jakob. Causal network reconstruction from time series: From theoretical
    assumptions to practical estimation, Chaos 28, 075310 (2018);
    doi: 10.1063/1.5025050
"""
logistic3(;u₀ = rand(3), r₁ = 4, r₂ = 4, r₃ = 4,
    σx = 0.05, σy = 0.05, σz = 0.05) = logistic3(u₀, r₁, r₂, r₃, σx, σy, σz)

function eom_logistic4(u, p, t)
    r₁, r₂, r₃, r₄, c₁₂, c₂₃, c₃₄  = (p...,)
    y₁, y₂, y₃, y₄ = (u...,)

    dy₁ = y₁*(r₁ - r₁*y₁)
    dy₂ = y₂*(r₂ - c₁₂*y₁ - r₂*y₂)
    dy₃ = y₃*(r₃ - c₂₃*y₂ - r₃*y₃)
    dy₄ = y₄*(r₄ - c₃₄*y₃ - r₄*y₄)
    return SVector{4}(dy₁, dy₂, dy₃, dy₄)
end

function logistic4(u₀, r₁, r₂, r₃, r₄, c₁₂, c₂₃, c₃₄)
    @warn "`logistic4` is deprecated in CausalityTools v2. "*
    "Use `system(Logistic4Chain())` instead, which returns a "*
    "`DiscreteDynamicalSystem` that can be iterated."
    p = @LArray [r₁, r₂, r₃, r₄, c₁₂, c₂₃, c₃₄] (:r₁, :r₂, :r₃, :r₄, :c₁₂, :c₂₃, :c₃₄)
    DiscreteDynamicalSystem(eom_logistic4, u₀, p)
end

"""
    logistic4(;u₀ = rand(4), r₁ = 3.9, r₂ = 3.6, r₃ = 3.6, r₄ = 3.8,
        c₁₂ = 0.4, c₂₃ = 0.4, c₃₄ = 0.35) → DiscreteDynamicalSystem

Initialise a system of a transitive chain of four unidirectionally coupled
logistic maps, where ``y_1 \\to y_2 \\to y_3 \\to y_4`` [1]. Default
parameters are as in [1].

*Note: With the default parameters which are as in [1], for some initial conditions,
this system wanders off to ``\\pm \\infty`` for some of the variables. Make sure that
you have a good realisation before using the orbit for anything.*

## Equations of motion

```math
\\begin{aligned}
y_1(t+1) &= y_1(t)(r_1 - r_1 y_1) \\\\
y_2(t+1) &= y_2(t)(r_2 - c_{12} y_1 - r_2 y_2) \\\\
y_3(t+1) &= y_3(t)(r_3 - c_{23} y_2 - r_3 y_3) \\\\
y_4(t+1) &= y_4(t)(r_4 - c_{34} y_3 - r_4 y_4)
\\end{aligned}
```

## References

1. Ye, Hao, et al. "Distinguishing time-delayed causal interactions using
    convergent cross mapping." Scientific reports 5 (2015): 14750
"""
logistic4(;u₀ = rand(4),
            r₁ = 3.9, r₂ = 3.6, r₃ = 3.6, r₄ = 3.8,
            c₁₂ = 0.4, c₂₃ = 0.4, c₃₄ = 0.35) =
    logistic4(u₀, r₁, r₂, r₃, r₄, c₁₂, c₂₃, c₃₄)


export nonlinear3d

"""
    eom_nonlinear3d(u₀, a₁, a₂, a₃,  b₁, b₂, b₃,
        c₁₂, c₂₃, c₁₃, σ₁, σ₂, σ₃) → DiscreteDynamicalSystem

Equations of motion for a 3d nonlinear system with nonlinear couplings
``x_1 \\to x_2``, ``x_2 \\to x_3`` and ``x_1 \\to x_3``. Modified from [1].

## Equations of motion

The equations of motion are

```math
\\begin{aligned}
x_1(t+1) &= a_1 x_1 (1-x_1(t))^2  e^{-x_2(t)^2} + 0.4 \\xi_{1}(t) \\\\
x_2(t+1) &= a_1 x_2 (1-x_2(t))^2  e^{-x_2(t)^2} + 0.4 \\xi_{2}(t) + b x_1 x_2 \\\\
x_3(t+1) &= a_3 x_3 (1-x_3(t))^2  e^{-x_3(t)^2} + 0.4 \\xi_{3}(t) + c x_{2}(t) + d x_{1}(t)^2.
\\end{aligned}
```

## References

1. Gourévitch, B., Le Bouquin-Jeannès, R., & Faucon, G. (2006). Linear and nonlinear
    causality between signals: methods, examples and neurophysiological
    applications. Biological Cybernetics, 95(4), 349–369.
"""
function eom_nonlinear3d(x, p, n)
    x₁, x₂, x₃ = (x...,)
    a₁, a₂, a₃, b₁, b₂, b₃, c₁₂, c₂₃, c₁₃, σ₁, σ₂, σ₃ = (p...,)
    ξ₁ = rand(Normal(0, σ₁))
    ξ₂ = rand(Normal(0, σ₂))
    ξ₃ = rand(Normal(0, σ₃))

    dx₁ = a₁*x₁*(1-x₁)^2 * exp(-x₁^2) + b₁*ξ₁
    dx₂ = a₂*x₂*(1-x₂)^2 * exp(-x₂^2) + b₂*ξ₂ + c₁₂*x₁*x₂
    dx₃ = a₃*x₃*(1-x₃)^2 * exp(-x₃^2) + b₃*ξ₃ + c₂₃*x₂ + c₁₃*x₁^2

    return SVector{3}(dx₁, dx₂, dx₃)
end

function nonlinear3d(u₀, a₁, a₂, a₃,  b₁, b₂, b₃, c₁₂, c₂₃, c₁₃, σ₁, σ₂, σ₃)
    @warn "`nonlinear3d` is deprecated in CausalityTools v2. "*
    "Use `system(Nonlinear3())` instead, which returns a "*
    "`DiscreteDynamicalSystem` that can be iterated."
    p = @LArray [a₁, a₂, a₃,  b₁, b₂, b₃, c₁₂, c₂₃, c₁₃, σ₁, σ₂, σ₃] (:a₁, :a₂, :a₃,  :b₁, :b₂, :b₃, :c₁₂, :c₂₃, :c₁₃, :σ₁, :σ₂, :σ₃)
    s = DiscreteDynamicalSystem(eom_nonlinear3d, u₀, p)
    return s
end

"""
    nonlinear3d(;u₀ = rand(3),
        σ₁ = 1.0, σ₂ = 1.0, σ₃ = 1.0,
        a₁ = 3.4, a₂ = 3.4, a₃ = 3.4,
        b₁ = 0.4, b₂ = 0.4, b₃ = 0.4,
        c₁₂ = 0.5, c₂₃ = 0.3, c₁₃ = 0.5) → DiscreteDynamicalSystem

A 3d nonlinear system with nonlinear couplings ``x_1 \\to x_2``,
``x_2 \\to x_3`` and ``x_1 \\to x_3``. Modified from [1].

## Equations of motion

The equations of motion are

```math
\\begin{aligned}
x_1(t+1) &= a_1 x_1 (1-x_1(t))^2  e^{-x_2(t)^2} + 0.4 \\xi_{1}(t) \\\\
x_2(t+1) &= a_1 x_2 (1-x_2(t))^2  e^{-x_2(t)^2} + 0.4 \\xi_{2}(t) + b x_1 x_2 \\\\
x_3(t+1) &= a_3 x_3 (1-x_3(t))^2  e^{-x_3(t)^2} + 0.4 \\xi_{3}(t) + c x_{2}(t) + d x_{1}(t)^2.
\\end{aligned}
```

## References

1. Gourévitch, B., Le Bouquin-Jeannès, R., & Faucon, G. (2006). Linear and nonlinear
    causality between signals: methods, examples and neurophysiological
    applications. Biological Cybernetics, 95(4), 349–369.
"""
nonlinear3d(;u₀ = rand(3),
        σ₁ = 1.0, σ₂ = 1.0, σ₃ = 1.0,
        a₁ = 3.4, a₂ = 3.4, a₃ = 3.4,
        b₁ = 0.4, b₂ = 0.4, b₃ = 0.4,
        c₁₂ = 0.5, c₂₃ = 0.3, c₁₃ = 0.5) =
    nonlinear3d(u₀, a₁, a₂, a₃,  b₁, b₂, b₃, c₁₂, c₂₃, c₁₃, σ₁, σ₂, σ₃)


function eom_nontrivial_pegiun(u, p, n)
    n = n + 10
    O = zeros(Float64, n + 3, 2)
    x, y = (u...,)
    p₁, p₂, p₃, p₄, p₅, p₆, σ₁, σ₂ = (p...,)

    # Propagate initial condition to the three first time steps.
    for i = 1:3
        O[i, 1] = x
        O[i, 2] = y
    end
    for i = 4:n
        y1 = O[i-1, 2]
        x2 = O[i-2, 1]
        y3 = O[i-3, 2]

        ξ₁ = rand(Normal(0, σ₁))
        ξ₂ = rand(Normal(0, σ₂))
        ynew = p₁*y1 + ξ₁
        xnew = p₂ + p₃*x2 + (p₄ - p₅*y3)/(1 + exp(-p₆*y3)) + ξ₂
        O[i, 1] = xnew
        O[i, 2] = ynew
    end
    O = O[10+3:end-10, :]
    O[:, 1] .= O[:, 1] .- mean(O[:, 1])
    O[:, 2] .= O[:, 2] .- mean(O[:, 2])
    O[:, 1] .= O[:, 1] ./ std(O[:, 1])
    O[:, 2] .= O[:, 2] ./ std(O[:, 2])
    return StateSpaceSet(O)
end

function nontrivial_pegiun(u₀, p₁, p₂, p₃, p₄, p₅, p₆, σ₁, σ₂, n::Int)
    @warn "`nontrivial_pegiun` is deprecated in CausalityTools v2. "*
    "Use `system(Peguin2())` instead, which returns a "*
    "`DiscreteDynamicalSystem` that can be iterated."
    p = @LArray [p₁, p₂, p₃, p₄, p₅, p₆, σ₁, σ₂] (:p₁, :p₂, :p₃, :p₄, :p₅, :p₆, :σ₁, :σ₂)
    eom_nontrivial_pegiun(u₀, p, n)
end


"""
    nontrivial_pegiun(;u₀ = rand(2), σ₁ = 0.1, σ₂ = 0.1,
        p₁ = 0.7, p₂ = 0.1, p₃ = 0.4, p₄ = 2.4, p₅ = 0.9, p₆ = 4, n = 100) → StateSpaceSet

A 2D discrete autoregressive system with nonlinear, nontrivial coupling from [1] .
This system is from [1](https://www.amse-aixmarseille.fr/sites/default/files/_dt/greqam/99a42.pdf), and
was also studied in [2](https://www.sciencedirect.com/science/article/pii/S0165027002003679).
The version implemented here allows for tweaking the parameters of the equations.
The difference equations are

```math
\\begin{aligned}
x(t+1) &= p_2 + p_3 x(t-2) + c_{yx}\\dfrac{p_4 - p_5 y(t-3)}{1 + e^{-p_6 y(t-3)}} + \\xi_1(t) \\
y(t+1) &= p_1 y(t) + \\xi_2(t).
\\end{aligned}
```
Here, ``\\xi_{1,2}(t)`` are two independent normally distributed noise processes
with zero mean and standard deviations ``\\sigma_1`` and ``\\sigma_2``. The
``\\xi_{1,2}(t)`` terms represent dynamical noise.

# References

[1] Péguin-Feissolle, A., & Teräsvirta, T. (1999). A General Framework for
Testing the Granger Noncausaality Hypothesis. Universites d’Aix-Marseille II
et III. [https://www.amse-aixmarseille.fr/sites/default/files/_dt/greqam/99a42.pdf](https://www.amse-aixmarseille.fr/sites/default/files/_dt/greqam/99a42.pdf)

[2] Chávez, M., Martinerie, J., & Le Van Quyen, M. (2003). Statistical
assessment of nonlinear causality: application to epileptic EEG signals.
Journal of Neuroscience Methods, 124(2), 113–128.
doi:10.1016/s0165-0270(02)00367-9
[https://www.sciencedirect.com/science/article/pii/S0165027002003679](https://www.sciencedirect.com/science/article/pii/S0165027002003679)
"""
function nontrivial_pegiun(;u₀ = rand(2), σ₁ = 0.1, σ₂ = 0.1,
        p₁ = 0.7, p₂ = 0.1, p₃ = 0.4, p₄ = 2.4, p₅ = 0.9, p₆ = 4, n = 100)
    eom_nontrivial_pegiun(u₀, [p₁, p₂, p₃, p₄, p₅, p₆, σ₁, σ₂], n)
end

function eom_ulam(dx, x, p, t)
    ε = p[:ε]
    f = x -> 2 - x^2
    dx[1] = f(ε*x[length(dx)] + (1-ε)*x[1])
    for i in 2:length(dx)
        dx[i] = f(ε*x[i-1] + (1-ε)*x[i])
    end
end

"""
    ulam(D::Int = 10; u₀ = rand(D), ε::Real = 0.10) → DiscreteDynamicalSystem

A lattice of `D` unidirectionally coupled ulam maps [Schreiber2000](@cite) defined as

```math
x^{m}_{t+1} = f(\\epsilon x^{m-1}_{t} + (1 - \\epsilon) x_{t}^{m}),
```

where ``m = 1, 2, \\ldots, D`` and ``f(x) = 2 - x^2``. In this system, information transfer
happens only in the direction of increasing ``m``.

"""
function ulam(D::Int = 10; u₀ = rand(D), ε::Real = 0.10)
    @warn "`ulam` is deprecated in CausalityTools v2. "*
    "Use `system(UlamLattice())` instead, which returns a "*
    "`DiscreteDynamicalSystem` that can be iterated."
    p = LVector(ε = ε)

    DiscreteDynamicalSystem(eom_ulam, u₀, p)
end

function eom_var1(x, p, n)
    σ₁, σ₂, σ₃ = p[1], p[2], p[3]
    x₁, x₂, x₃ = x[1], x[2], x[3]
    θ = rand(Normal(0, σ₁))
    η = rand(Normal(0, σ₂))
    ϵ = rand(Normal(0, σ₃))

    dx₁ = θ
    dx₂ = x₁ * η
    dx₃ = 0.5*x₃ * x₂ + ϵ
    return SVector{3}(dx₁, dx₂, dx₃)
end

function var1(u₀, σ₁, σ₂, σ₃)
    @warn "`var1` is deprecated in CausalityTools v2. "*
    "Use `system(Var1())` instead, which returns a "*
    "`DiscreteDynamicalSystem` that can be iterated."
    p = LVector(ε = ε)
    p = @LArray [σ₁, σ₂, σ₃] (:σ₁, :σ₂, :σ₃)

    DiscreteDynamicalSystem(eom_var1, u₀, p)
end

"""
    var1(x, p, n) → DiscreteDynamicalSystem

Initialise a discrete vector autoregressive system where X₁ → X₂ → X₃.
"""
var1(;u₀ = rand(3), σ₁ = 1.0, σ₂ = 0.2, σ₃ = 0.3) = var1(u₀, σ₁, σ₂, σ₃)

function eom_verdes(u, p, t)
    x, y, z = (u...,)
    ωy, ωz, σx, σy, σz = (p...,)

    ηx = σx == 0 ? 0 : rand(Normal(0, σx))
    ηy = σy == 0 ? 0 : rand(Normal(0, σy))
    ηz = σz == 0 ? 0 : rand(Normal(0, σz))

    dx = y*(18y - 27y^2 + 10)/2 + z*(1-z) + ηx
    dy = (1 - cos((2*pi/ωy) * t))/2 + ηy
    dz = (1 - sin((2*pi/ωz) * t))/2 + ηz
    return SVector{3}(dx, dy, dz)
end

function verdes(u₀, ωy, ωz, σx, σy, σz)
    @warn "`verdes` is deprecated in CausalityTools v2. "*
    "Use `system(Verdes())` instead, which returns a "*
    "`DiscreteDynamicalSystem` that can be iterated."
    p = LVector(ωy = ωy, ωz = ωz, σx = σx, σy = σy, σz = σz)

    DiscreteDynamicalSystem(eom_verdes, u₀, p)
end

"""
    verdes(;u₀ = rand(3), ωy = 315, ωz = 80,
        σx = 0.0, σy = 0.0, σz = 0.0) → DiscreteDynamicalSystem

Intitialise a 3D system where the response X is a highly nonlinear combination
of Y and Z [Verdes2005](@cite). The forcings Y and Z involve sines and cosines, and
have different periods, which controlled by `ωy` and `ωz`.

The equations of motion are

```math
\\begin{aligned}
x(t+1) &= \\dfrac{y(t)(18y(t) - 27y(t)^2 + 10)}{2} + z(t)(1-z(t)) + ηx \\
y(t+1) &= \\dfrac{(1 - \\dfrac{\\cos(2\\pi)}{\\omega y}t)}{2} + ηy \\
z(t+1) &= \\dfrac{(1 - \\dfrac{\\sin(2\\pi)}{\\omega z}t)}{2} + ηz
\\end{aligned}
```
where ηx, ηy, ηz is gaussian noise with mean 0 and standard deviation `σx`, `σy`
and `σz`.
"""
verdes(;u₀ = rand(3),
    ωy = 315, ωz = 80,
    σx = 0.01, σy = 0.01, σz = 0.01) =
    verdes(u₀, ωy, ωz, σx, σy, σz)

import Distributions: Distribution, Uniform, Normal

"""
    noise_uu(n::Int, lo = - 1, hi = 1)

Generate a signal consisting of `n` steps of uncorrelated uniform noise from
a uniform distribution on `[lo, hi]`.
"""
function noise_uu(n::Int; lo = - 1, hi = 1)
    @warn "`noise_uu` is deprecated in CausalityTools v2. "
    u = Uniform(-lo, hi)
    rand(u, n)
end


"""
    noise_ug(n::Int; μ = 0, σ = 1)

Generate a signal consisting of `n` steps of uncorrelated Gaussian noise from
a normal distribution with mean `μ` and standard deviation `σ`.
"""
function noise_ug(n::Int; μ = 0, σ = 1)
    @warn "`noise_ug` is deprecated in CausalityTools v2. "
    d = Normal(μ, σ)
    rand(d, n)
end

"""
    noise_brownian(n::Int; lo = - 1, hi = 1)
    noise_brownian(d::Distribution, n::Int)

Generate a signal consisting of `n` steps of Brownian noise, generated as
the zero-mean and unit standard deviation normalised cumulative sum of noise
generated from a uniform distribution on `[lo, hi]`. Optionally, a distribution
`d` from which to sample can be provided.

## Examples

```julia
# Based on uncorrelated uniform noise
noise_brownian(100)
noise_brownian(100, lo = -2, hi = 2)
noise_brownian(Uniform(-3, 3), 100)

# Based on uncorrelated Gaussian noise
μ, σ = 0, 2
noise_brownian(Normal(μ, σ), 100)
```
"""
function noise_brownian(n::Int; lo = - 1, hi = 1)
    @warn "`noise_brownian` is deprecated in CausalityTools v2."

    u = Uniform(lo, hi)
    xs = cumsum(rand(u, n))
    (xs .- mean(xs)) ./ std(xs)
end

function noise_brownian(d::Distribution, n::Int)
    @warn "`noise_brownian` is deprecated in CausalityTools v2."
    xs = cumsum(rand(d, n))
    (xs .- mean(xs)) ./ (std(xs))
end

export noise_uu, noise_ug, noise_brownian
