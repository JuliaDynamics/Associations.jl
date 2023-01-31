export independence, cit
export ConditionalIndependenceTest

abstract type IndependenceTest end
"""
    ConditionalIndependenceTest <: IndependenceTest

The supertype for all conditional independence tests, which are:

- [`LocalPermutation`](@ref).
"""
abstract type ConditionalIndependenceTest <: IndependenceTest end

"""
    independence(test, x, y, [z])

Perform the given `test` of independence between `x` and `y` using the provided `test`.
If `z` is given, compute the conditional independence of `x` and `y` given `z`.

This function is just a generic implementation of a one-sided hypothesis test, where the
null hypothesis is that `x` and `y` are independent (given `z`, if provided).

## Supported tests

The null hypothesis is specified by `test`, which is a [`IndependenceTest`](@ref).

- [`LocalPermutation`](@ref)
- [`SurrogateCIT`](@ref). This is essentially a convenience wrapper that performs
    standard surrogate testing.
"""
function independence(test, args...; kwargs...)
    error("No concrete implementation for $(typeof(test)) test yet")
end

include("local_permutation/LocalPermutation.jl")
include("surrogate/SurrogateTest.jl")
