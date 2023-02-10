export independence, cit
export IndependenceTest


"""
    IndependenceTest <: IndependenceTest

The supertype for all independence tests.
"""
abstract type IndependenceTest end

"""
    independence(test::IndependenceTest, x, y, [z]) → summary

Perform the given [`IndependenceTest`](@ref) `test` on data `x`, `y` and `z`.
If only `x` and `y` are given, `test` must provide a bivariate association measure.
If `z` is given too, then `test` must provide a conditional association measure.

Returns a test `summary`, whose type depends on `test`.

## Compatible independence tests

- [`SurrogateTest`](@ref).
- [`LocalPermutationTest`](@ref).
- [`JointDistanceDistributionTest`](@ref).
"""
function independence(test, args...; kwargs...)
    error("No concrete implementation for $(typeof(test)) test yet")
end

include("parametric/parametric.jl")
include("local_permutation/LocalPermutationTest.jl")
include("surrogate/SurrogateTest.jl")
# TODO: rename/find suitable generic name before including
# include("correlation/correlation.jl). Perhaps `ParametricTest`?
