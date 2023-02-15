using Random
using TimeseriesSurrogates
export SurrogateTest
export SurrogateTestResult

"""
    SurrogateTest <: IndependenceTest
    SurrogateTest(measure, [est];
        nshuffles::Int = 100,
        surrogate = RandomShuffle(),
        rng = Random.default_rng(),
    )

A generic (conditional) independence test for assessing whether two variables `X` and `Y`
are independendent, potentially conditioned on a third variable `Z`, based on
surrogate data. Used with [`independence`](@ref).

## Description

This is a generic one-sided hypothesis test that checks whether `x` and `y`
are independent (given `z`, if provided) based on resampling from a null distribution
assumed to represent independence between the variables. The null distribution is generated
by repeatedly shuffling the input data in some way that is intended
to break any dependence between the input variables.

There are different ways of shuffling, dictated by `surrogate`, each representing a
distinct null hypothesis. For each shuffle, the provided `measure` is computed (using `est`,
if relevant). This procedure is repeated `nshuffles` times, and a test summary is returned.
The shuffled variable is always the first variable (`X`). Exceptions are:

- If [`TransferEntropy`](@ref) measure such as [`TEShannon`](@ref),
    then the source variable is always shuffled, and the target and conditional
    variable are left unshuffled.

## Compatible measures

| Measure                               | Pairwise | Conditional | Requires `est` |
| ------------------------------------- | :------: | :---------: | :------------: |
| [`PearsonCorrelation`](@ref)          |    ✓    |     ✖      |       No       |
| [`DistanceCorrelation`](@ref)         |    ✓    |     ✖      |       No       |
| [`SMeasure`](@ref)                    |    ✓    |     ✖      |       No       |
| [`HMeasure`](@ref)                    |    ✓    |     ✖      |       No       |
| [`MMeasure`](@ref)                    |    ✓    |     ✖      |       No       |
| [`LMeasure`](@ref)                    |    ✓    |     ✖      |       No       |
| [`PairwiseAsymmetricInference`](@ref) |    ✓    |     ✖      |      Yes       |
| [`ConvergentCrossMapping`](@ref)      |    ✓    |     ✖      |      Yes       |
| [`MIShannon`](@ref)                   |    ✓    |     ✖      |      Yes       |
| [`MIRenyiJizba`](@ref)                |    ✓    |     ✖      |      Yes       |
| [`MIRenyiSarbu`](@ref)                |    ✓    |     ✖      |      Yes       |
| [`MITsallisMartin`](@ref)             |    ✓    |     ✖      |      Yes       |
| [`MITsallisFuruichi`](@ref)           |    ✓    |     ✖      |      Yes       |
| [`PartialCorrelation`](@ref)          |    ✖    |     ✓      |      Yes       |
| [`CMIShannon`](@ref)                  |    ✖    |     ✓      |      Yes       |
| [`CMIRenyiJizba`](@ref)               |    ✖    |     ✓      |      Yes       |
| [`TransferEntropy`](@ref)             |    ✓    |     ✓      |      Yes       |

## Examples

- [Pairwise test, `DistanceCorrelation`](@ref examples_surrogatetest_distancecorrelation).
- [Pairwise test, `TEShannon`](@ref examples_surrogatetest_teshannon).
- [Conditional test, `PartialCorrelation`](@ref examples_surrogatetest_partialcorrelation).
- [Pairwise test, `MIShannon`, categorical](@ref examples_surrogatetest_mishannon_categorical).
- [Conditional test, `CMIShannon`, categorical](@ref examples_surrogatetest_cmishannon_categorical).
"""
struct SurrogateTest{M, E, R, S} <: IndependenceTest
    measure::M
    est::E
    rng::R
    surrogate::S
    nshuffles::Int

    function SurrogateTest(measure::M, est::E = nothing;
        rng::R = Random.default_rng(),
        surrogate::S = RandomShuffle(),
        nshuffles::Int = 100,
        ) where {M, E, R, S}
        new{M, E, R, S}(measure, est, rng, surrogate, nshuffles)
    end
end


Base.show(io::IO, test::SurrogateTest) = print(io,
    """
    `SurrogateTest` independence test.
    -------------------------------------
    measure:    $(test.measure)
    estimator:  $(test.est)
    rng:        $(test.rng)
    # shuffles: $(test.nshuffles)
    surrogate:  $(test.surrogate)
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
    nshuffles::Int
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
        Ensemble quantiles ($(test.nshuffles) permutations):
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
    (; measure, est, rng, surrogate, nshuffles) = test
    X, Y, Z = Dataset(x), Dataset(y), Dataset(z)
    @assert length(X) == length(Y) == length(Z)
    N = length(x)
    Î = estimate(measure,est, X, Y, Z)
    s = surrogenerator(x, surrogate, rng)
    Îs = zeros(nshuffles)
    for b in 1:nshuffles
        Îs[b] = estimate(measure, est, s(), Y, Z)
    end
    p = count(Î .<= Îs) / nshuffles

    return SurrogateTestResult(Î, Îs, p, nshuffles)
end

function independence(test::SurrogateTest, x, y)
    (; measure, est, rng, surrogate, nshuffles) = test
    X, Y = Dataset(x), Dataset(y)
    @assert length(X) == length(Y)
    N = length(x)
    Î = estimate(measure,est, X, Y)
    sx = surrogenerator(x, surrogate, rng)
    Îs = zeros(nshuffles)
    for b in 1:nshuffles
        Îs[b] = estimate(measure, est, sx(), y)
    end
    p = count(Î .<= Îs) / nshuffles

    return SurrogateTestResult(Î, Îs, p, nshuffles)
end

# Concrete implementations
include("contingency.jl")
include("transferentropy.jl")
include("hlms_measure.jl")
include("crossmapping.jl")
include("mutualinfo.jl")
include("condmutualinfo.jl")
