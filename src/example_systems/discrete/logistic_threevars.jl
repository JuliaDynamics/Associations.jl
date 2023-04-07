using Random: default_rng
using DynamicalSystemsBase: DiscreteDynamicalSystem
using Graphs: add_edge!, edges
using Graphs.SimpleGraphs: SimpleDiGraph
using StaticArrays: SVector
using SimpleWeightedGraphs: SimpleWeightedGraph

export Uncoupled
export SingleUnidir
export CommonCause, CommonCauseSingle, CommonCauseMutual
export MutualSingle, MutualDouble, MutualTriple
export DirectedChain

Base.@kwdef struct Uncoupled{V, C, Σ, RNG} <: DiscreteDefinition
    xi::V = rand(3)
    c::C = 0.3
    σ::Σ = 0.2
    f::Function = (x) -> 4.0*(x - x^2)
    g::Function = (a, b, c, σ, rng) -> (b + c*(a + σ * rand(rng)) ) / (1 + c*(1+σ))
    h::Function = (a, b) -> (a/2 + b/2) % 1.0
    rng::RNG = default_rng()
end

function system(definition::Uncoupled)
    return DiscreteDynamicalSystem(eom, definition.xi, definition)
end

function eom(u, p::Uncoupled, t)
    (; xi, c, σ, f, g, h, rng) = p
    x, y, z = u
    dx = f(x)
    dy = f(y)
    dz = f(z)
    return SVector{3}(dx, dy, dz)
end

function SimpleWeightedDiGraph(sys::Uncoupled)
    g = SimpleWeightedDiGraph(3)
    if sys.c > 0
        add_edge!(g, 3, 1,  sys.c)
        add_edge!(g, 3, 2,  sys.c)
    end
    return g
end

function SimpleDiGraph(sys::Uncoupled)
    g = SimpleDiGraph(3)
    if sys.c > 0
        add_edge!(g, 3, 1)
        add_edge!(g, 3, 2)
    end
    return g
end

Base.@kwdef struct SingleUnidir{V, C, Σ, RNG} <: DiscreteDefinition
    xi::V = rand(3)
    c::C = 0.3
    σ::Σ = 0.2
    f::Function = (x) -> 4.0*(x - x^2)
    g::Function = (a, b, c, σ, rng) -> (b + c*(a + σ * rand(rng)) ) / (1 + c*(1+σ))
    h::Function = (a, b) -> (a/2 + b/2) % 1.0
    rng::RNG = default_rng()
end

function system(definition::SingleUnidir)
    return DiscreteDynamicalSystem(eom, definition.xi, definition)
end

function eom(u, p::SingleUnidir, t)
    (; xi, c, σ, f, g, h, rng) = p
    x, y, z = u
    dx = f(x)
    dy = g(x, y, c, σ, rng) |> f
    dz = f(z)
    return SVector{3}(dx, dy, dz)
end

function SimpleWeightedDiGraph(sys::SingleUnidir)
    g = SimpleWeightedDiGraph(3)
    if sys.c > 0
        add_edge!(g, 1, 2, sys.c)
    end
    return g
end

function SimpleDiGraph(sys::SingleUnidir)
    g = SimpleDiGraph(3)
    if sys.c > 0
        add_edge!(g, 1, 2)
    end
    return g
end

Base.@kwdef struct CommonCause{V, C, Σ, RNG} <: DiscreteDefinition
    xi::V = rand(3)
    c::C = 0.3
    σ::Σ = 0.2
    f::Function = (x) -> 4.0*(x - x^2)
    g::Function = (a, b, c, σ, rng) -> (b + c*(a + σ * rand(rng)) ) / (1 + c*(1+σ))
    h::Function = (a, b) -> (a/2 + b/2) % 1.0
    rng::RNG = default_rng()
end

function system(definition::CommonCause)
    return DiscreteDynamicalSystem(eom, definition.xi, definition)
end

function eom(u, p::CommonCause, t)
    (; xi, c, σ, f, g, h, rng) = p
    x, y, z = u
    dx = g(z, x, c, σ, rng) |> f
    dy = g(z, y, c, σ, rng) |> f
    dz = f(z)
    return SVector{3}(dx, dy, dz)
end

function SimpleWeightedDiGraph(sys::CommonCause)
    g = SimpleWeightedDiGraph(3)
    if sys.c > 0
        add_edge!(g, 3, 1,  sys.c)
        add_edge!(g, 3, 2,  sys.c)
    end
    return g
end

function SimpleDiGraph(sys::CommonCause)
    g = SimpleDiGraph(3)
    if sys.c > 0
        add_edge!(g, 3, 1)
        add_edge!(g, 3, 2)
    end
    return g
end


Base.@kwdef struct CommonCauseSingle{V, C, Σ, RNG} <: DiscreteDefinition
    xi::V = rand(3)
    c::C = 0.3
    σ::Σ = 0.2
    f::Function = (x) -> 4.0*(x - x^2)
    g::Function = (a, b, c, σ, rng) -> (b + c*(a + σ * rand(rng)) ) / (1 + c*(1+σ))
    h::Function = (a, b) -> (a/2 + b/2) % 1.0
    rng::RNG = default_rng()
end

function system(definition::CommonCauseSingle)
    return DiscreteDynamicalSystem(eom, definition.xi, definition)
end

function eom(u, p::CommonCauseSingle, t)
    (; xi, c, σ, f, g, h, rng) = p
    x, y, z = u
    xy = g(x, y, c, σ, rng) |> f
    yx = g(y, x, c, σ, rng) |> f
    zy = g(z, y, c, σ, rng) |> f

    dx = g(z, x, c, σ, rng) |> f
    dy = h(xy, zy)
    dz = f(z)
    return SVector{3}(dx, dy, dz)
end

function SimpleWeightedDiGraph(sys::CommonCauseSingle)
    g = SimpleWeightedDiGraph(3)
    if sys.c > 0
        add_edge!(g, 1, 2,  sys.c)
        add_edge!(g, 3, 2,  sys.c)
        add_edge!(g, 3, 1,  sys.c)
    end
    return g
end

function SimpleDiGraph(sys::CommonCauseSingle)
    g = SimpleDiGraph(3)
    if sys.c > 0
        add_edge!(g, 1, 2)
        add_edge!(g, 3, 2)
        add_edge!(g, 3, 1)
    end
    return g
end


Base.@kwdef struct CommonCauseMutual{V, C, Σ, RNG} <: DiscreteDefinition
    xi::V = rand(3)
    c::C = 0.3
    σ::Σ = 0.2
    f::Function = (x) -> 4.0*(x - x^2)
    g::Function = (a, b, c, σ, rng) -> (b + c*(a + σ * rand(rng)) ) / (1 + c*(1+σ))
    h::Function = (a, b) -> (a/2 + b/2) % 1.0
    rng::RNG = default_rng()
end

function system(definition::CommonCauseMutual)
    return DiscreteDynamicalSystem(eom, definition.xi, definition)
end

function eom(u, p::CommonCauseMutual, t)
    (; xi, c, σ, f, g, h, rng) = p
    x, y, z = u
    xy = g(x, y, c, σ, rng) |> f
    yx = g(y, x, c, σ, rng) |> f
    zx = g(z, x, c, σ, rng) |> f
    zy = g(z, y, c, σ, rng) |> f
    dx = h(zx, yx)
    dy = h(zy, xy)
    dz = f(z)
    return SVector{3}(dx, dy, dz)
end

function SimpleWeightedDiGraph(sys::CommonCauseMutual)
    g = SimpleWeightedDiGraph(3)
    if sys.c > 0
        add_edge!(g, 1, 2,  sys.c)
        add_edge!(g, 2, 1,  sys.c)
        add_edge!(g, 3, 1,  sys.c)
        add_edge!(g, 3, 2,  sys.c)
    end
    return g
end

function SimpleDiGraph(sys::CommonCauseMutual)
    g = SimpleDiGraph(3)
    if sys.c > 0
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 1)
        add_edge!(g, 3, 1)
        add_edge!(g, 3, 2)
    end
    return g
end

################################################################

Base.@kwdef struct MutualSingle{V, C, Σ, RNG} <: DiscreteDefinition
    xi::V = rand(3)
    c::C = 0.3
    σ::Σ = 0.2
    f::Function = (x) -> 4.0*(x - x^2)
    g::Function = (a, b, c, σ, rng) -> (b + c*(a + σ * rand(rng)) ) / (1 + c*(1+σ))
    h::Function = (a, b) -> (a/2 + b/2) % 1.0
    rng::RNG = default_rng()
end

function system(definition::MutualSingle)
    return DiscreteDynamicalSystem(eom, definition.xi, definition)
end

function eom(u, p::MutualSingle, t)
    (; xi, c, σ, f, g, h, rng) = p
    x, y, z = u
    dx = g(y, x, c, σ, rng) |> f
    dy = g(x, y, c, σ, rng) |> f
    dz = f(z)
    return SVector{3}(dx, dy, dz)
end

function SimpleWeightedDiGraph(sys::MutualSingle)
    g = SimpleWeightedDiGraph(3)
    if sys.c > 0
        add_edge!(g, 1, 2,  sys.c)
        add_edge!(g, 2, 1,  sys.c)
    end
    return g
end

function SimpleDiGraph(sys::MutualSingle)
    g = SimpleDiGraph(3)
    if sys.c > 0
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 1)
    end
    return g
end

Base.@kwdef struct MutualDouble{V, C, Σ, RNG} <: DiscreteDefinition
    xi::V = rand(3)
    c::C = 0.3
    σ::Σ = 0.2
    f::Function = (x) -> 4.0*(x - x^2)
    g::Function = (a, b, c, σ, rng) -> (b + c*(a + σ * rand(rng)) ) / (1 + c*(1+σ))
    h::Function = (a, b) -> (a/2 + b/2) % 1.0
    rng::RNG = default_rng()
end

function system(definition::MutualDouble)
    return DiscreteDynamicalSystem(eom, definition.xi, definition)
end

function eom(u, p::MutualDouble, t)
    (; xi, c, σ, f, g, h, rng) = p
    x, y, z = u
    xy = g(x, y, c, σ, rng) |> f
    yx = g(y, x, c, σ, rng) |> f
    zy = g(z, y, c, σ, rng) |> f
    yz = g(y, z, c, σ, rng) |> f
    dx = yx
    dy = h(zy, xy)
    dz = yz
    return SVector{3}(dx, dy, dz)
end

function SimpleWeightedDiGraph(sys::MutualDouble)
    g = SimpleWeightedDiGraph(3)
    if sys.c > 0
        add_edge!(g, 1, 2,  sys.c)
        add_edge!(g, 2, 1,  sys.c)
        add_edge!(g, 2, 3,  sys.c)
        add_edge!(g, 3, 2,  sys.c)
    end
    return g
end

function SimpleDiGraph(sys::MutualDouble)
    g = SimpleDiGraph(3)
    if sys.c > 0
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 1)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 2)
    end
    return g
end

Base.@kwdef struct MutualTriple{V, C, Σ, RNG} <: DiscreteDefinition
    xi::V = rand(3)
    c::C = 0.3
    σ::Σ = 0.2
    f::Function = (x) -> 4.0*(x - x^2)
    g::Function = (a, b, c, σ, rng) -> (b + c*(a + σ * rand(rng)) ) / (1 + c*(1+σ))
    h::Function = (a, b) -> (a/2 + b/2) % 1.0
    rng::RNG = default_rng()
end

function system(definition::MutualTriple)
    return DiscreteDynamicalSystem(eom, definition.xi, definition)
end

function eom(u, p::MutualTriple, t)
    (; xi, c, σ, f, g, h, rng) = p
    x, y, z = u
    xy = g(x, y, c, σ, rng) |> f
    yx = g(y, x, c, σ, rng) |> f
    zy = g(z, y, c, σ, rng) |> f
    yz = g(y, z, c, σ, rng) |> f
    zx = g(z, x, c, σ, rng) |> f
    xz = g(x, z, c, σ, rng) |> f
    dx = h(zx, yx)
    dy = h(zy, xy)
    dz = h(xz, yz)
    return SVector{3}(dx, dy, dz)
end

function SimpleWeightedDiGraph(sys::MutualTriple)
    g = SimpleWeightedDiGraph(3)
    if sys.c > 0
        add_edge!(g, 1, 2,  sys.c)
        add_edge!(g, 2, 1,  sys.c)
        add_edge!(g, 1, 3,  sys.c)
        add_edge!(g, 3, 1,  sys.c)
        add_edge!(g, 2, 3,  sys.c)
        add_edge!(g, 3, 2,  sys.c)
    end
    return g
end

function SimpleDiGraph(sys::MutualTriple)
    g = SimpleDiGraph(3)
    if sys.c > 0
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 1)
        add_edge!(g, 1, 3)
        add_edge!(g, 3, 1)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 2)
    end
    return g
end

Base.@kwdef struct DirectedChain{V, C, Σ, RNG} <: DiscreteDefinition
    xi::V = rand(3)
    c::C = 0.3
    σ::Σ = 0.2
    f::Function = (x) -> 4.0*(x - x^2)
    g::Function = (a, b, c, σ, rng) -> (b + c*(a + σ * rand(rng)) ) / (1 + c*(1+σ))
    h::Function = (a, b) -> (a/2 + b/2) % 1.0
    rng::RNG = default_rng()
end

function system(definition::DirectedChain)
    return DiscreteDynamicalSystem(eom, definition.xi, definition)
end

function eom(u, p::DirectedChain, t)
    (; xi, c, σ, f, g, h, rng) = p
    x, y, z = u
    dx = g(z, x, c, σ, rng) |> f
    dy = g(x, y, c, σ, rng) |> f
    dz = g(y, z, c, σ, rng) |> f
    return SVector{3}(dx, dy, dz)
end

function SimpleWeightedDiGraph(sys::DirectedChain)
    g = SimpleWeightedDiGraph(3)
    if sys.c > 0
        add_edge!(g, 1, 2,  sys.c)
        add_edge!(g, 2, 3,  sys.c)
        add_edge!(g, 3, 1,  sys.c)
    end
    return g
end

function SimpleDiGraph(sys::DirectedChain)
    g = SimpleDiGraph(3)
    if sys.c > 0
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 1)
    end
    return g
end
