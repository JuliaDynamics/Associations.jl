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

A surrogate-data based generic (conditional) independence test for assessing whether the 
association between variables `X` and `Y` are independent, potentially conditioned on a 
third variable `Z`.

## Compatible estimators and measures 

- Compatible with [`AssociationMeasure`](@ref)s that measure some sort of pairwise or conditional association.

!!! note 
    You must yourself determine whether using a particular measure is *meaningful*, and what it
    *means*.

!!! note 
    If used with a [`TransferEntropy`](@ref) measure such as [`TEShannon`](@ref),
    then the source variable is always shuffled, and the target and conditional
    variable are left unshuffled.

## Usage 

- Use with [`independence`](@ref) to perform a surrogate test with input data. This will
    return a [`SurrogateAssociationTestResult`](@ref).

## Description

This is a generic one-sided hypothesis test that checks whether `x` and `y`
are independent (given `z`, if provided) based on resampling from a null distribution
assumed to represent independence between the variables. The null distribution is generated
by repeatedly shuffling the input data in some way that is intended
to break any dependence between the input variables.

The test first estimates the desired statistic using `est_or_measure` on the input data. 
Then, the first input variable is shuffled `nshuffled` times according to the given 
`surrogate` method (each type of `surrogate` represents a distinct null hypothesis).
For each shuffle, `est_or_measure` is recomputed and the results are stored. 

## Examples

- [Example 1](@ref example_SurrogateAssociationTest_SMeasure):
     [`SMeasure`](@ref) test for pairwise independence.
- [Example 2](@ref example_SurrogateAssociationTest_DistanceCorrelation): 
    [`DistanceCorrelation`](@ref) test for pairwise independence.
- [Example 3](@ref example_SurrogateAssociationTest_PartialCorrelation):
    [`PartialCorrelation`](@ref) test for conditional independence.
- [Example 4](@ref example_SurrogateAssociationTest_MIShannon_categorical):
    [`MIShannon`](@ref) test for pairwise independence on categorical data.
- [Example 5](@ref example_SurrogateAssociationTest_CMIShannon_categorical):
    [`CMIShannon`](@ref) test for conditional independence on categorical data.  
- [Example 6](@ref example_independence_MCR): [`MCR`](@ref) test for 
    pairwise and conditional independence.  
- [Example 7](@ref example_SurrogateAssociationTest_ChatterjeeCorrelation).
    [`ChatterjeeCorrelation`](@ref) test for pairwise independence.
- [Example 8](@ref example_SurrogateAssociationTest_AzadkiaChatterjeeCoefficient).
    [`AzadkiaChatterjeeCoefficient`](@ref) test for pairwise and conditional independence.
"""
struct SurrogateAssociationTest{E, R, S} <: IndependenceTest{E}
    est_or_measure::E
    rng::R
    surrogate::S
    nshuffles::Int
    show_progress::Bool
end
function SurrogateAssociationTest(
        est_or_measure::E;
        rng::R = Random.default_rng(),
        surrogate::S = RandomShuffle(),
        nshuffles::Int = 100, 
        show_progress = false
        ) where {E, R, S}
    SurrogateAssociationTest{E, R, S}(est_or_measure, rng, surrogate, nshuffles, show_progress)
end


Base.show(io::IO, test::SurrogateAssociationTest) = print(io,
    """
    `SurrogateAssociationTest` independence test.
    -------------------------------------
    estimator/measure: $(test.est_or_measure)
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
    (; est_or_measure, rng, surrogate, nshuffles, show_progress) = test
    verify_number_of_inputs_vars(est_or_measure, 1+length(args))
    SSSets = map(w -> StateSpaceSet(w), args)
    estimation = x -> association(est_or_measure, x, SSSets...)
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
    return SurrogateAssociationTestResult(1+length(args), Î, Îs, p, nshuffles)
end

# Concrete implementations
include("transferentropy.jl")
include("crossmapping.jl")
include("hlms_measure.jl")
include("chatterjee_correlation.jl")
include("azadkia_chatterjee_coefficient.jl")

# Input checks
function SurrogateAssociationTest(measure::T) where T <: MultivariateInformationMeasure
    str = "`SurrogateAssociationTest` can't be constructed using the information measure `$T` definition directly. " * 
        "Give a valid estimator as the first argument instead and give the " * 
        "definition to the estimator, e.g. " * 
        "FPVP(CMIShannon())"
    throw(ArgumentError(str))
end