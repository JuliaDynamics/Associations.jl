export independence
export ConditionalIndependenceTest

abstract type IndependenceTest end
"""
    ConditionalIndependenceTest <: IndependenceTest

The supertype for all conditional independence tests.

Currently, the following concrete implementations exist:

- [`LocalPermutation`](@ref).
"""
abstract type ConditionalIndependenceTest <: IndependenceTest end

"""
    independence(test, x, y, z)

Compute the conditional independence of `x` and `y` given `z` using
the provided [`ConditionalIndependenceTest`](@ref).

## Implemented independence tests

- [`LocalPermutation`](@ref)
"""
function independence(test, args...; kwargs...)
    error("No concrete implementation for $(typeof(test)) test yet")
end

include("LocalPermutation.jl")
