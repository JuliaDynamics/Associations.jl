using TimeseriesSurrogates

"""
    GlobalPermutation <: ConditionalIndependenceTest


The `GlobalPermutation` test is a generic test that uses the surrogate data technique to
assess whether two variables `X` and `Y` are conditionally independendent given a third
variable `Z`.

The shuffled variable is `X` (the first input to [`independence`](@ref)). Since most
surrogate data are defined for timeseries only, so the `GlobalPermutation` test also
will mostly just work when `X` is a univariate timeseries.

See also: [`LocalPermutation`](@ref),
[TimeseriesSurrogates.jl](https://github.com/JuliaDynamics/TimeseriesSurrogates.jl)
"""
Base.@kwdef struct GlobalPermutation{D, M, E, R, S} <: ConditionalIndependenceTest
    definition::D = CMI4HShannon()
    measure::M = CMIShannon(base = 2)
    est::EST = VejmelkaPalus(k = 5)
    rng::R = MersenneTwister(12345678)
    surrogate::S = RandomShuffle()
    nsurr::Int = 100
end

function independence(test::GlobalPermutation, x, y, z)
    (; definition, measure, est, rng, surrogate, nsurr) = test
    X, Y, Z = Dataset(x), Dataset(y), Dataset(z)
    @assert length(X) == length(Y) == length(Z)
    N = length(x)
    Î = estimate(measure,est, X, Y, Z)
    s = surrogenerator(x, surrogate, rng)
    Îs = zeros(nsurr)
    for b in 1:nsurr
        x̂ = s()
        Îs[b] = estimate(definition, measure, est, x̂, Y, Z)
    end
    p = count(Î .<= Îs) / nsurr

    return GlobalPermutationTest(Î, Îs, p, nsurr)
end
