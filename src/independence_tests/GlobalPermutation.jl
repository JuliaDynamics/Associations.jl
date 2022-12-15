"""
    GlobalPermutation <: ConditionalIndependenceTest

The `GlobalPermutation` conditional independence test uses surrogate
data to shuffle.

A related test is [`LocalPermutation`](@ref), but for that test, shuffled
data preserve the local

See also:
[TimeseriesSurrogates.jl](https://github.com/JuliaDynamics/TimeseriesSurrogates.jl)
"""
Base.@kwdef struct GlobalPermutation{M <: CMI, E <: Entropy, EST, R, S} <: ConditionalIndependenceTest
    measure::M = CMI()
    e::E = Shannon(; base = 2)
    est::EST = VejmelkaPalus(k = 5)
    rng::R = Random.default_rng()
    surrogate::S = 5
    nsurr::Int = 100
end
