using LabelledArrays: @LArray
using StaticArrays: SVector
using DynamicalSystemsBase: DiscreteDynamicalSystem
using Random

export ChaoticNoisyLinear2

"""
    ChaoticNoisyLinear2 <: DiscreteDefinition
    ChaoticNoisyLinear2(; xi = [0.1, 0.2], c = 0.5,
        nx = Normal(0, 0.05), ny = Normal(0, 0.05),
        rng = Random.default_rng())

A bivariate system of two chaotic maps that are linearly coupled from `x → y`
with coupling strength `c`.

## Definition

```math
\\begin{align*}
x(t+1) = 3.4 x(t) (1 - x(t)^2) e^{-x(t)^2} + 0.8x(t-1) + \\xi_x \\\\
y(t+1) = 3.4 y(t) (1 - y(t)^2) e^{-y(t)^2} + 0.8y(t-1) + \\xi_y + c x(t-2)
\\end{align*}
```
Process noise ``\\xi_x`` and ``\\xi_y``
is drawn at each iteration from `nx` and `ny`.
"""
struct ChaoticNoisyLinear2{P, V, NX, NY, C, R} <: LaggedDiscreteDefinition{P}
    past_states::P
    xi::V
    nx::NX
    ny::NY
    c::C
    rng::R

    function ChaoticNoisyLinear2(; xi::V = [0.1, 0.2], c::C = 0.5,
            nx::NX = Normal(0, 0.05),
            ny::NY = Normal(0, 0.05),
            rng::R = Random.default_rng()) where {V, C, NX, NY, R}
        T = eltype(1.0)
        mx = MVector{2, T}(repeat([xi[1]], 2))
        my = MVector{2, T}(repeat([xi[2]], 2))
        past_states = SVector{2, MVector{2, T}}(mx, my)
        P = typeof(past_states)
        return new{P, V, NX, NY, C, R}(past_states, xi, nx, ny, c, rng)
    end
end

function system(definition::ChaoticNoisyLinear2)
    return DiscreteDynamicalSystem(eom_linearmap2, definition.xi, definition)
end

function eom_linearmap2(u, p::ChaoticNoisyLinear2, t)
    (; past_states, xi, nx, ny, c, rng) = p
    # `u` is simply ignored here, because the state is stored in the memory vectors
    mx, my = past_states
    x₁, x₂ = mx[1], mx[2]
    y₁, y₂ = my[1], my[2]
    dx = 3.4 * x₁ * (1 - x₁^2) * exp(-x₁^2) + 0.8*x₂ + rand(rng, nx)
    dy = 3.4 * y₁ * (1 - y₁^2) * exp(-y₁^2) + 0.5*y₂ + c*x₂ + rand(rng, ny)
    new_state = SVector{2}(dx, dy)
    update_states!(p, new_state) # Update memory vectors
    return new_state
end
