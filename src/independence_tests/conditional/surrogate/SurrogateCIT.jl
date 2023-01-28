using TimeseriesSurrogates
export SurrogateCIT
export SurrogateCITResult

"""
    SurrogateCIT <: IndependenceTest
    SurrogateCIT(;
        measure = CMIShannon(),
        est = FPVP(k = 5),
        nsurr::Int = 100,
        surrogate = RandomShuffle(),
        rng = Random.MersenneTwister(1234),
    )

The `SurrogateCIT` test is a generic conditional independence test (CIT) for assessing
whether two variables `X` and `Y` are conditionally independendent given a third variable
`Z`. It uses
[surrogate time series](https://github.com/JuliaDynamics/TimeseriesSurrogates.jl) to
generate the null distribution.

## Null hypotheses

- If used with a [`ConditionalMutualInformation`](@ref) measure such as
    [`CMIShannon`](@ref), then the shuffled variable is the first input `X`.
- If used with a [`TransferEntropy`](@ref) measure such as
    [`TEShannon`](@ref), then the shuffled variable is the first input `X`, i.e.
    the source variable.


## Example usage

Here, we'll test if the transfer entropy from `x` to `y` is significant. If we can reject
the null hypothesis that `x` and `y` are independent, then we take that as evidence
that `x` influences `y`.

```julia
using CausalityTools
# If `x` and `y` were dependent given `z`, then we'd expect to be able to reject the
# null hypothesis that x ⫫ y | z.
sys = logistic2_unidir(c_xy = 0.5)
npts = 1000
x, y = columns(trajectory(sys, npts, Ttr = npts*10))
test = SurrogateCIT(TEShannon(), FPVP(k = 10))
independence(test, x, y)

# If `x` and `y` are *independent*.
```
"""
struct SurrogateCIT{M, E, R, S} <: IndependenceTest
    measure::M
    est::E
    rng::R
    surrogate::S
    nsurr::Int

    function SurrogateCIT(measure::M = CMIShannon(base = 2), est::EST = FPVP(k = 5);
        rng::R = MersenneTwister(1234),
        surrogate::S = RandomShuffle(),
        nsurr::Int = 100,
        ) where {M, EST, R, S}
        new{M, EST, R, S}(measure, est, rng, surrogate, nsurr)
    end
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
        `SurrogateCIT` independence test result
        -----------------------------------------------------------------------------------
        H₀: "The first two variables are independent (given the 3rd variable, if relevant)"
        Hₐ: "The first two variables are dependent (given the 3rd variable, if relevant)"
        -----------------------------------------------------------------------------------
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

# Generic dispatch for any three-argument conditional independence measure where the
# third argument is to be conditioned on. This works naturally with e.g.
# conditional mutual information.
function independence(test::SurrogateCIT, x, y, z)
    (; measure, est, rng, surrogate, nsurr) = test
    X, Y, Z = Dataset(x), Dataset(y), Dataset(z)
    @assert length(X) == length(Y) == length(Z)
    N = length(x)
    Î = estimate(measure,est, X, Y, Z)
    s = surrogenerator(x, surrogate, rng)
    Îs = zeros(nsurr)
    for b in 1:nsurr
        Îs[b] = estimate(measure, est, s(), Y, Z)
    end
    p = count(Î .<= Îs) / nsurr

    return SurrogateCITResult(Î, Îs, p, nsurr)
end

include("transferentropy.jl")
