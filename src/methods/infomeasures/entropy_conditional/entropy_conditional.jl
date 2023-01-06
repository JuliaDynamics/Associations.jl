export entropy_conditional
export ConditionalEntropy
export ConditionalEntropyDefinition

"""
The supertype for all conditional entropies.
"""
abstract type ConditionalEntropy <: InformationMeasure end

"""
The supertype for all conditional entropy definitions.
"""
abstract type ConditionalEntropyDefinition <: Definition end

# Measures
include("CEShannon.jl")
include("CETsallisFuruichi.jl")
include("CETsallisAbe.jl")

"""
    entropy_conditional(measure::ConditionalEntropy, c::ContingencyMatrix{T, 2}) where T

Compute the given conditional entropy `measure` from the pre-computed
[`ContingencyMatrix`](@ref) `c`.

## Supported measures

- [`CETsallis`](@ref). Tsallis conditional entropy.

"""
entropy_conditional(measure::ConditionalEntropy, args...; kwargs...) =
    estimate(measure, args...; kwargs...)
