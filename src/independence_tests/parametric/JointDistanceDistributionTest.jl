export JointDistanceDistributionTest
export JDDTestResult

using Statistics: mean, std
using Distributions: TDist
using Distributions: cdf

"""
    JointDistanceDistributionTest <: IndependenceTest
    JointDistanceDistributionTest(measure::JointDistanceDistribution; rng = Random.default_rng())

An independence test for two variables based on the [`JointDistanceDistribution`](@ref)
(Amigó & Hirata, 2018)[^Amigo2018].

When used with [`independence`](@ref), a [`JDDTestResult`](@ref) is returned.

## Description

The joint distance distribution (labelled `Δ` in their paper) is used by Amigó & Hirata (2018)
to detect directional couplings of the form ``X \\to Y`` or ``Y \\to X``. `JointDistanceDistributionTest`
formulates their method as an independence test.

Formally, we test the hypothesis ``H_0`` (the variables are independent) against ``H_1``
(there is directional coupling between the variables). To do so, we use a right-sided/upper-tailed t-test
to check mean of `Δ` is skewed towards positive value, i.e.

* ``H_0 := \\mu(\\Delta) = 0``
* ``H_1 := \\mu(\\Delta) > 0``.

When used with [`independence`](@ref), a `JDDTestResult` is returned, which contains
the joint distance distribution and a p-value. If you only need `Δ`, use [`jdd`](@ref) directly.

## Examples

[This example](@ref quickstart_jddtest) shows how the `JointDistanceDistributionTest`
can be used in practice.

[^Amigo2018]:
    Amigó, José M., and Yoshito Hirata. "Detecting directional couplings from multivariate
    flows by the joint distance distribution." Chaos: An Interdisciplinary Journal of
    Nonlinear Science 28.7 (2018): 075302.
"""
struct JointDistanceDistributionTest{M <: JointDistanceDistribution, R} <: IndependenceTest{M}
    measure::M
    rng::R

    function JointDistanceDistributionTest(
            measure::M = JointDistanceDistribution();
            rng::R = Random.default_rng()) where {M, R}
        new{M, R}(measure, rng)
    end
end

# Performs a one-sided t-test to see if the joint distance distribution is skewed above
# `measure.μ`, which is the hypothetical mean of the joint distance
# distribution under the null (defaults to `0.0`).
"""
    JDDTestResult(Δjdd, hypothetical_μ, pvalue)

Holds the results of [`JointDistanceDistributionTest`](@ref). `Δjdd` is the
`Δ`-distribution, `hypothetical_μ` is the hypothetical mean of the `Δ`-distribution
under the null, and `pvalue` is the p-value for the one-sided t-test.

## Implements

- **`pvalue`**. Returns a p-value for the test.
- **`point_estimate`**. Returns the mean of the estimated `Δ` distribution.
"""
struct JDDTestResult{V, T, P} <: IndependenceTestResult
    n_vars::Int # 2 vars = pairwise, 3 vars = conditional
    Δjdd::V
    hypothetical_μ::T
    pvalue::P
end

pvalue(x::JDDTestResult) = x.pvalue
point_estimate(x::JDDTestResult) = mean(x.Δjdd)

function Base.show(io::IO, r::JDDTestResult)
    # TODO: make a function to do this (a pairwise and a conditional version), so this
    # isn't repeated everywhere.
    α005 = r.pvalue < 0.05 ?
        "α = 0.05:  ✓ Evidence favors dependence" :
        "α = 0.05:  ✖ Independence cannot be rejected"
    α001 = r.pvalue < 0.01 ?
        "α = 0.01:  ✓ Evidence favors dependence" :
        "α = 0.01:  ✖ Independence cannot be rejected"
    α0001 = r.pvalue < 0.001 ?
        "α = 0.001: ✓ Evidence favors dependence" :
        "α = 0.001: ✖ Independence cannot be rejected"

    print(io,
        """\
        `JointDistanceDistributionTest` independence test
        ----------------------------------------------------------------------------------
        H₀: μ(Δ) = 0 (the input variables are independent)
        H₁: μ(Δ) > 0 (there is directional dependence between the input variables)
        ----------------------------------------------------------------------------------
        Hypothetical μ(Δ): $(r.hypothetical_μ)
        Observed μ(Δ):     $(mean(r.Δjdd))
        p-value:           $(r.pvalue)
          $α005
          $α001
          $α0001\
        """
        )
end


function independence(test::JointDistanceDistributionTest, x, y)
    Δjdd = jdd(test.measure, x, y)

    # Right-sided t-test
    t = t_statistic(Δjdd, hypothetical_μ = test.measure.μ)
    degrees_of_freedom = length(Δjdd) - 1
    D = TDist(degrees_of_freedom)
    pval = 1 - cdf(D, t)

    return JDDTestResult(2, Δjdd, test.measure.μ, pval)
end
