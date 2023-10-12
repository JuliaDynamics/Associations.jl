export ConditionalEntropy

"""
    ConditionalEntropy <: MultivariateInformationMeasure

The supertype for all conditional entropy measures.

## Concrete subtypes

- [`ConditionalEntropyShannon`](@ref)
- [`ConditionalEntropyTsallisAbe`](@ref)
- [`ConditionalEntropyTsallisFuruichi`](@ref)
"""
abstract type ConditionalEntropy <: MultivariateInformationMeasure end

min_inputs_vars(::ConditionalEntropy) = 2
max_inputs_vars(::ConditionalEntropy) = 2

include("ConditionalEntropyShannon.jl")
include("ConditionalEntropyTsallisAbe.jl")
include("ConditionalEntropyTsallisFuruichi.jl")