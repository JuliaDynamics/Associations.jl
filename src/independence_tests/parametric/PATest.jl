export PATest
export PATestResult

"""
    PATest <: IndependenceTest
    PATest(measure::PA, est)

A test for directional (conditional) dependence based on the modified [`PA`](@ref) measure,
estimated using the given estimator `est`. Compatible estimators are listed in
the docstring for [`PA`](@ref).

When used with [`independence`](@ref), a [`PATestResult`](@ref) summary is returned.

!!! note
    This is an experimental test. It is part of an ongoing paper submission revision,
    but is provided here for convenience.
"""
struct PATest{M <: PA, E <: PA_ESTS} <: IndependenceTest
    measure::M
    est::E
end

"""
    PATestResult(n_vars, ΔA, ttest, pvalue)

Holds the result of a [`SurrogateTest`](@ref). `n_vars` is the number of variables
used for the test (2 for pairwise, 3 for conditional). `ΔA` is the distribution of
asymmetries, one for each `η`. `ttest` is a one-sample t-test, and `pvalue` is the
right-tailed p-value for the test.
"""
struct PATestResult{A, T, P} <: IndependenceTestResult
    n_vars::Int
    ΔA::A
    ttest::T
    pvalue::P
end

function independence(test::PATest, x, y)
    ΔA = estimate(test.measure, test.est, x, y)
    t = OneSampleTTest(ΔA, 0)
    p = pvalue(t, tail = :right)
    return PATestResult(2, ΔA, t, p)
end

function independence(test::PATest, x, y, z)
    ΔA = estimate(test.measure, test.est, x, y, z)
    t = OneSampleTTest(ΔA, 0)
    p = pvalue(t, tail = :right)
    return PATestResult(3, ΔA, t, p)
end


pvalue(r::PATestResult) = r.pvalue

function Base.show(io::IO, test::PATestResult)
    print(io,
        """\
        `PATest` independence test
        $(null_hypothesis_text(test))
        $(pvalue_text_summary(test))
        """
        )
end
