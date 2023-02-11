using LabelledArrays: @LArray
using StaticArrays: SVector
using DynamicalSystemsBase: DiscreteDynamicalSystem
using Random

export LinearMap2

struct LinearMap2{V, T, NX, NY, C, R} <: DiscreteDefinition
    xi::V
    mx::MVector{2, T} # holds past states of x1
    my::MVector{2, T} # holds past states of x2
    nx::NX
    ny::NY
    c::C
    rng::R

    function LinearMap2(; xi::V = [0.1, 0.2], c::C = 0.5,
            nx::NX = Normal(0, 0.05),
            ny::NY = Normal(0, 0.05),
            rng::R = Random.default_rng()) where {V, C, NX, NY, R}
        T = eltype(1.0)
        mx = MVector{2, T}(repeat([xi[1]], 2))
        my = MVector{2, T}(repeat([xi[2]], 2))
        return new{V, T, NX, NY, C, R}(xi, mx, my, nx, ny, c, rng)
    end
end

function system(definition::LinearMap2)
    return DiscreteDynamicalSystem(eom_linearmap2, definition.xi, definition)
end

function eom_linearmap2(u, p::LinearMap2, t)
    (; xi, mx, my, nx, ny, c, rng) = p
    # `u` is simply ignored here, because the state is stored in the memory vectors
    x₁, x₂ = mx[1], mx[2]
    y₁, y₂ = my[1], my[2]
    dx = 3.4 * x₁ * (1 - x₁^2) * exp(-x₁^2) + 0.8*x₂ + rand(rng, nx)
    dy = 3.4 * y₁ * (1 - y₁^2) * exp(-y₁^2) + 0.5*y₂ + c*x₂ + rand(rng, ny)
    new_state = SVector{2}(dx, dy)
    update_state!(p, new_state) # Update memory vectors
    return new_state
end

function update_state!(p::LinearMap2, xnew::SVector{2})
    p.mx[2] = p.mx[1]
    p.mx[1] = xnew[1]
    p.my[2] = p.my[1]
    p.my[1] = xnew[2]
end
