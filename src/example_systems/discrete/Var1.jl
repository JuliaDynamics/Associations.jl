using DynamicalSystemsBase: DiscreteDynamicalSystem
using Distributions: Normal
using Random

export Var1

"""
    Var1 <: DiscreteDefinition
    Var1(; xi = [0.5, 0.5, 0.5],
        a = 0.5, θ = Normal(0, 1.0), η = Normal(0, 0.2), ϵ = Normal(0, 0.3),
        rng = Random.default_rng())

A discrete vector autoregressive system where X₁ → X₂ → X₃.
"""
Base.@kwdef struct Var1{V, A, Σ, N, E, R} <: DiscreteDefinition
    xi::V = [0.5, 0.5, 0.5]
    a::A = 0.5
    θ::Σ = Normal(0, 1.0)
    η::N = Normal(0, 0.2)
    ϵ::E = Normal(0, 0.3)
    rng::R = Random.default_rng()
end

function system(definition::Var1)
    return DiscreteDynamicalSystem(eom_var1system, definition.xi, definition)
end

function eom_var1system(u, p::Var1, n)
    x₁, x₂, x₃ = u
    (; a, θ, η, ϵ, rng) = p
    dx₁ = rand(rng, θ)
    dx₂ = x₁ * rand(rng, η)
    dx₃ = a*x₃ * x₂ + rand(rng, ϵ)
    return SVector{3}(dx₁, dx₂, dx₃)
end
