using LabelledArrays: @LArray
using StaticArrays: SVector
using DynamicalSystemsBase: DiscreteDynamicalSystem
using DynamicalSystemsBase: trajectory
using Distributions: Normal
using Statistics: mean, std

export Peguin2

"""
    Peguin2 <: DiscreteDefinition
    Peguin2(; xi = [0.5, 0.4], σ₁ = 0.1, σ₂ = 0.1,
        p₁ = 0.7, p₂ = 0.1, p₃ = 0.4, p₄ = 2.4, p₅ = 0.9, p₆ = 4)

A 2D discrete autoregressive system with nonlinear, nontrivial coupling from [1] .
This system is from Péguin-Feissolle & Teräsvirta (1999)[^Péguin-Feissolle1999], and
was also studied in Chávez et al. (2003)[^Chávez2003].

## Description

The system is defined by the equations

```math
\\begin{align*}
x(t+1) &= p_2 + p_3 x(t-2) + c_{yx}\\dfrac{p_4 - p_5 y(t-3)}{1 + e^{-p_6 y(t-3)}} + \\xi_1(t) \\\\
y(t+1) &= p_1 y(t) + \\xi_2(t).
\\end{align*}
```

Here, ``\\xi_{1,2}(t)`` are two independent normally distributed noise processes
with zero mean and standard deviations ``\\sigma_1`` and ``\\sigma_2``. The
``\\xi_{1,2}(t)`` terms represent dynamical noise. The parameters of the original system
are here tweakable.

[^Péguin-Feissolle1999]:
    Péguin-Feissolle, A., & Teräsvirta, T. (1999). A General Framework for
    Testing the Granger Noncausaality Hypothesis. Universites d’Aix-Marseille II
    et III. [https://www.amse-aixmarseille.fr/sites/default/files/_dt/greqam/99a42.pdf](https://www.amse-aixmarseille.fr/sites/default/files/_dt/greqam/99a42.pdf)

[^Chávez2003]:
    Chávez, M., Martinerie, J., & Le Van Quyen, M. (2003). Statistical
    assessment of nonlinear causality: application to epileptic EEG signals.
    Journal of Neuroscience Methods, 124(2), 113–128.
    doi:10.1016/s0165-0270(02)00367-9
    [https://www.sciencedirect.com/science/article/pii/S0165027002003679](https://www.sciencedirect.com/science/article/pii/S0165027002003679)
"""
struct Peguin2{P,V,T,Nx,Ny,P1,P2,P3,P4,P5,P6,R} <: LaggedDiscreteDefinition{P}
     # `past_states[i]` := past states of the i-th variable of the system.
    past_states::P
    xi::V
    nx::Nx # a distribution to sample noise from for x
    ny::Ny # a distribution to sample noise from for y
    p₁::P1
    p₂::P2
    p₃::P3
    p₄::P4
    p₅::P5
    p₆::P6
    rng::R

    function Peguin2(;
            xi::V = [0.4, 0.5],
            nx::Nx = Normal(0, 0.1),
            ny::Ny = Normal(0, 0.1),
            p₁::P1 = 0.7,
            p₂::P2 = 0.1,
            p₃::P3 = 0.4,
            p₄::P4 = 2.4,
            p₅::P5 = 0.9,
            p₆::P6 = 4.0,
            rng::R = Random.default_rng()) where {V,Nx,Ny,P1,P2,P3,P4,P5,P6,R}
        T = eltype(1.0)

        mx = MVector{3, T}(repeat([xi[1]], 3))
        my = MVector{3, T}(repeat([xi[2]], 3))
        past_states = SVector{2, MVector{3, T}}(mx, my)
        P = typeof(past_states)
        return new{P, V,T,Nx,Ny,P1,P2,P3,P4,P5,P6,R}(
            past_states, xi, nx, ny, p₁, p₂, p₃, p₄, p₅, p₆, rng)
    end
end

function system(definition::Peguin2)
    return DiscreteDynamicalSystem(eom_peguin2, definition.xi, definition)
end

function eom_peguin2(u, p, t)
    (; xi, nx, ny, p₁, p₂, p₃, p₄, p₅, p₆, rng) = p
    # `u` is simply ignored here, because the state is stored in the memory vectors
    mx, my = p.past_states
    x₂ = mx[2]
    y₁, y₃ = my[1], my[3]
    dx = p₁*y₁ + rand(rng, nx)
    dy = p₂ + p₃*x₂ + (p₄ - p₅*y₃)/(1 + exp(-p₆*y₃)) + rand(rng, ny)
    new_state = SVector{2}(dx, dy)
    @inbounds update_states!(p, new_state)
    return new_state
end
