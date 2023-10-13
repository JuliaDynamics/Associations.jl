using Random
using TimeseriesSurrogates
import ProgressMeter
export SurrogateAssociationTest
export SurrogateAssociationTestResult
import Statistics: quantile

"""
    SurrogateAssociationTest <: IndependenceTest
    SurrogateAssociationTest(est_or_measure;
        nshuffles::Int = 100,
        surrogate = RandomShuffle(),
        rng = Random.default_rng(),
        show_progress = false,
    )

A generic (conditional) independence test for assessing whether two variables `X` and `Y`
are independendent, potentially conditioned on a third variable `Z`, based on
surrogate data.

When used with [`independence`](@ref), a [`SurrogateAssociationTestResult`](@ref) is returned.

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

## Compatible estimators/measures

Some measures can be used directly as the first input, since they don't require any 
estimator, e.g. one can construct `SurrogateAssociationTest(PearsonCorrelation())`.

| Measure                       | Pairwise | Conditional |
| ----------------------------- | :------: | :---------: |
| [`PearsonCorrelation`](@ref)  |    ✓    |     ✖      |
| [`DistanceCorrelation`](@ref) |    ✓    |     ✓      |
| [`SMeasure`](@ref)            |    ✓    |     ✖      |
| [`HMeasure`](@ref)            |    ✓    |     ✖      |
| [`MMeasure`](@ref)            |    ✓    |     ✖      |
| [`LMeasure`](@ref)            |    ✓    |     ✖      |

Moreover, any valid estimator of the following measures may be used:

| Measure                               | Pairwise | Conditional |
| ------------------------------------- | :------: | :---------: |
| [`PairwiseAsymmetricInference`](@ref) |    ✓    |     ✖      |
| [`ConvergentCrossMapping`](@ref)      |    ✓    |     ✖      |
| [`MIShannon`](@ref)                   |    ✓    |     ✖      |
| [`MIRenyiJizba`](@ref)                |    ✓    |     ✖      |
| [`MIRenyiSarbu`](@ref)                |    ✓    |     ✖      |
| [`MITsallisMartin`](@ref)             |    ✓    |     ✖      |
| [`MITsallisFuruichi`](@ref)           |    ✓    |     ✖      |
| [`PartialCorrelation`](@ref)          |    ✖    |     ✓      |
| [`CMIShannon`](@ref)                  |    ✖    |     ✓      |
| [`CMIRenyiJizba`](@ref)               |    ✖    |     ✓      |
| [`TEShannon`](@ref)                   |    ✓    |     ✓      |
| [`TERenyiJizba`](@ref)                |    ✓    |     ✓      |
| [`PMI`](@ref)                         |    ✖    |     ✓      |



## Examples

- [Pairwise test, `DistanceCorrelation`](@ref examples_SurrogateAssociationTest_distancecorrelation).
- [Pairwise test, `TEShannon`](@ref examples_SurrogateAssociationTest_teshannon).
- [Conditional test, `PartialCorrelation`](@ref examples_SurrogateAssociationTest_partialcorrelation).
- [Pairwise test, `MIShannon`, categorical](@ref examples_SurrogateAssociationTest_mishannon_categorical).
- [Conditional test, `CMIShannon`, categorical](@ref examples_SurrogateAssociationTest_cmishannon_categorical).
"""
struct SurrogateAssociationTest{M, E, R, S} <: IndependenceTest{M}
    est_or_measure::M
    est::E
    rng::R
    surrogate::S
    nshuffles::Int
    show_progress::Bool
end
function SurrogateAssociationTest(est_or_measure::M, est::E = nothing;
    rng::R = Random.default_rng(),
    surrogate::S = RandomShuffle(),
    nshuffles::Int = 100, show_progress = false
    ) where {M, E, R, S}
    SurrogateAssociationTest{M, E, R, S}(est_or_measure, est, rng, surrogate, nshuffles, show_progress)
end


Base.show(io::IO, test::SurrogateAssociationTest) = print(io,
    """
    `SurrogateAssociationTest` independence test.
    -------------------------------------
    estimator/measure: $(test.est_or_measure)
    estimator:         $(test.est)
    rng:               $(test.rng)
    # shuffles:        $(test.nshuffles)
    surrogate:         $(test.surrogate)
    """
)

"""
    SurrogateAssociationTestResult(m, m_surr, pvalue)

Holds the result of a [`SurrogateAssociationTest`](@ref). `m` is the measure computed on
the original data. `m_surr` is a vector of the measure computed on permuted data, where
`m_surr[i]` is the measure compute on the `i`-th permutation. `pvalue` is the one-sided
`p`-value for the test.
"""
struct SurrogateAssociationTestResult{M, MS, P} <: IndependenceTestResult
    n_vars::Int # 2 vars = pairwise, 3 vars = conditional
    m::M
    m_surr::MS
    pvalue::P
    nshuffles::Int
end
pvalue(r::SurrogateAssociationTestResult) = r.pvalue
quantile(r::SurrogateAssociationTestResult, q) = quantile(r.m_surr, q)

function Base.show(io::IO, test::SurrogateAssociationTestResult)
    print(io,
        """\
        `SurrogateAssociationTest` independence test
        $(null_hypothesis_text(test))
        $(quantiles_text(test))
        $(pvalue_text_summary(test))
        """
        )
end

# Generic dispatch for any three-argument conditional independence measure where the
# third argument is to be conditioned on. This works naturally with e.g.
# conditional mutual information.
function independence(test::SurrogateAssociationTest, x, args...)
    # Setup (`args...` is either `y` or `y, z`)
    (; est_or_measure, est, rng, surrogate, nshuffles, show_progress) = test
    verify_number_of_inputs_vars(est_or_measure, 1+length(args))
    SSSets = map(w -> StateSpaceSet(w), args)
    estimation = x -> estimate(est_or_measure, x, SSSets...)
    progress = ProgressMeter.Progress(nshuffles;
        desc="SurrogateAssociationTest:", enabled=show_progress
    )

    # Estimate
    Î = estimation(StateSpaceSet(x))
    s = surrogenerator(x, surrogate, rng)
    Îs = zeros(nshuffles)
    for b in 1:nshuffles
        Îs[b] = estimation(s())
        ProgressMeter.next!(progress)
    end
    p = count(Î .<= Îs) / nshuffles
    return SurrogateAssociationTestResult(3, Î, Îs, p, nshuffles)
end

# Concrete implementations
# include("contingency.jl")
# include("transferentropy.jl")
# include("hlms_measure.jl")
# include("crossmapping.jl")
# include("pmi.jl")
