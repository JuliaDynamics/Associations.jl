using TimeseriesSurrogates
export SurrogateTest
export SurrogateTestResult

"""
    SurrogateTest <: IndependenceTest
    SurrogateTest(;
        measure = CMIShannon(),
        est = FPVP(k = 5),
        nsurr::Int = 100,
        surrogate = RandomShuffle(),
        rng = Random.MersenneTwister(1234),
    )

The `SurrogateTest` test is a generic conditional independence test (CIT) for assessing
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
test = SurrogateTest(TEShannon(), FPVP(k = 10))
independence(test, x, y)

# If `x` and `y` are *independent*.
```
"""
struct SurrogateTest{M, E, R, S} <: IndependenceTest
    measure::M
    est::E
    rng::R
    surrogate::S
    nsurr::Int

    function SurrogateTest(measure::M, est::E = nothing;
        rng::R = MersenneTwister(1234),
        surrogate::S = RandomShuffle(),
        nsurr::Int = 100,
        ) where {M, E, R, S}
        new{M, E, R, S}(measure, est, rng, surrogate, nsurr)
    end
end


Base.show(io::IO, test::SurrogateTest) = print(io,
    """
    `SurrogateTest` independence test.
    -------------------------------------
    measure:        $(test.measure)
    estimator:      $(test.est)
    rng:            $(test.rng)
    # permutations: $(test.nsurr)
    surrogate:      $(test.surrogate)
    """
)

"""
    SurrogateTestResult(M, Msurr, pvalue)

Holds the result of a [`SurrogateTestResult`](@ref). `M` is the measure computed on
the original data. `Msurr` is a vector of the measure computed on permuted data, where
Msurr[i] corresponds to the `i`-th permutation. `pvalue` is the `p`-value for the test.
"""
struct SurrogateTestResult{M, MS, P}
    M::M
    Msurr::MS
    pvalue::P
    nsurr::Int
end
pvalue(r::SurrogateTestResult) = r.pvalue
quantile(r::SurrogateTestResult, q) = quantile(r.Msurr, q)

function Base.show(io::IO, test::SurrogateTestResult)
    α005 = pvalue(test) < 0.05 ?
        "α = 0.05:  ✓ Evidence favors dependence" :
        "α = 0.05:  ✖ Independence cannot be rejected"
    α001 = pvalue(test) < 0.01 ?
        "α = 0.01:  ✓ Evidence favors dependence" :
        "α = 0.01:  ✖ Independence cannot be rejected"
    α0001 = pvalue(test) < 0.001 ?
        "α = 0.001: ✓ Evidence favors dependence" :
        "α = 0.001: ✖ Independence cannot be rejected"

    print(io,
        """\
        `SurrogateTest` independence test result
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
function independence(test::SurrogateTest, x, y, z)
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

    return SurrogateTestResult(Î, Îs, p, nsurr)
end

function independence(test::SurrogateTest, x, y)
    (; measure, est, rng, surrogate, nsurr) = test
    X, Y = Dataset(x), Dataset(y)
    @assert length(X) == length(Y)
    N = length(x)
    Î = estimate(measure,est, X, Y)
    sx = surrogenerator(x, surrogate, rng)
    sy = surrogenerator(y, surrogate, rng)
    Îs = zeros(nsurr)
    for b in 1:nsurr
        Îs[b] = estimate(measure, est, sx(), sy())
    end
    p = count(Î .<= Îs) / nsurr

    return SurrogateTestResult(Î, Îs, p, nsurr)
end

include("contingency.jl")
include("transferentropy.jl")
