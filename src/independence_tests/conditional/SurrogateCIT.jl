using TimeseriesSurrogates
export SurrogateCIT
export SurrogateCITResult

"""
    SurrogateCIT <: ConditionalIndependenceTest
    SurrogateCIT(;
        measure = CMIShannon(),
        est = FrenzelPompeVelmejkaPalus(k = 5),
        nsurr::Int = 100,
        surrogate = RandomShuffle(),
        rng = Random.MersenneTwister(1234),
    )

The `SurrogateCIT` test is a generic conditional independence test (CIT) that uses the
surrogate data technique to assess whether two variables `X` and `Y` are conditionally
independendent given a third variable `Z`.

The shuffled variable is `X` (the first input to [`independence`](@ref)). Since most
surrogate data are defined for timeseries only, so the `SurrogateCIT` test also
will mostly just work when `X` is a univariate timeseries.

See also: [`LocalPermutation`](@ref),
[TimeseriesSurrogates.jl](https://github.com/JuliaDynamics/TimeseriesSurrogates.jl)
"""
Base.@kwdef struct SurrogateCIT{M, E, R, S} <: ConditionalIndependenceTest
    measure::M = CMIShannon(base = 2)
    est::E = FrenzelPompeVelmejkaPalus(k = 5)
    rng::R = MersenneTwister(12345678)
    surrogate::S = RandomShuffle()
    nsurr::Int = 100
end


Base.show(io::IO, test::SurrogateCIT) = print(io,
    """
    `SurrogateCIT` independence test.
    -------------------------------------
    measure:        $(test.measure)
    estimator:      $(test.est)
    rng:            $(test.rng)
    # permutations: $(test.nsurr)
    surrogate:      $(test.surrogate)
    """
)

"""
    SurrogateCITResult(M, Msurr, pvalue)

Holds the result of a [`SurrogateCITResult`](@ref). `M` is the measure computed on
the original data. `Msurr` is a vector of the measure computed on permuted data, where
Msurr[i] corresponds to the `i`-th permutation. `pvalue` is the `p`-value for the test.
"""
struct SurrogateCITResult{M, MS, P}
    M::M
    Msurr::MS
    pvalue::P
    nsurr::Int
end
pvalue(r::SurrogateCITResult) = r.pvalue
quantile(r::SurrogateCITResult, q) = quantile(r.Msurr, q)

function Base.show(io::IO, test::SurrogateCITResult)
    α005 = pvalue(test) < 0.05 ?
        "H₀ rejected at α = 0.05:  Yes ✅" :
        "H₀ rejected at α = 0.05:  No  ❌"
    α001 = pvalue(test) < 0.01 ?
        "H₀ rejected at α = 0.01:  Yes ✅" :
        "H₀ rejected at α = 0.01:  No  ❌"
    α0001 = pvalue(test) < 0.001 ?
        "H₀ rejected at α = 0.001: Yes ✅" :
        "H₀ rejected at α = 0.001: No  ❌"

    print(io,
        """\
        `SurrogateCIT` independence test
        ----------------------------------------------------------------------------------
        H₀: "The first two variables are conditionally independent given the third"
        ----------------------------------------------------------------------------------
        Estimated: $(test.M)
        Ensemble quantiles ($(test.nsurr) permutations):
          (99.9%): $(quantile(test.Msurr, 0.999))
          (99%):   $(quantile(test.Msurr, 0.99))
          (95%):   $(quantile(test.Msurr, 0.95))
        p-value:   $(test.pvalue)
          $α005
          $α001
          $α0001\
        """
        )
end

function conditional_independence(test::SurrogateCIT, x, y, z)
    (; measure, est, rng, surrogate, nsurr) = test
    X, Y, Z = Dataset(x), Dataset(y), Dataset(z)
    @assert length(X) == length(Y) == length(Z)
    N = length(x)
    Î = estimate(measure,est, X, Y, Z)
    s = surrogenerator(x, surrogate, rng)
    Îs = zeros(nsurr)
    for b in 1:nsurr
        x̂ = s()
        Îs[b] = estimate(measure, est, x̂, Y, Z)
    end
    p = count(Î .<= Îs) / nsurr

    return SurrogateCITResult(Î, Îs, p, nsurr)
end


# function conditional_independence(test::SurrogateCIT{M <: TransferEntropy}, x...)
#     (; measure, est, rng, surrogate, nsurr) = test

#     S, T, Tf, C = individual_marginals(measure.embedding, x...)
#     cmi = te_to_cmi(measure)
#     # TE(s -> t) := I(t⁺; s⁻ | t⁻, c⁻).
#     return condmutualinfo(cmi, est, Tf, S, Dataset(T, C))

#     T⁺, S⁻, T⁻ = Dataset(t⁺), Dataset(s⁻), Dataset(t⁻)
#     @assert length(T⁺) == length(S⁻) == length(T⁻)
#     N = length(x)
#     Î = transferentropy(measure.est, T⁺, Y, Z)
#     s = surrogenerator(x, surrogate, rng)
#     Îs = zeros(nsurr)
#     for b in 1:nsurr
#         x̂ = s()
#         Îs[b] = estimate(measure, est, x̂, Y, Z)
#     end
#     p = count(Î .<= Îs) / nsurr

#     return SurrogateCITResult(Î, Îs, p, nsurr)
# end
