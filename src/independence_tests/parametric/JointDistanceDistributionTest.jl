export JointDistanceDistributionTest
using Statistics: mean

"""
    JointDistanceDistributionTest <: IndependenceTest
    JointDistanceDistributionTest(measure::JointDistanceDistribution;
        rng = Random.default_rng(),
        n_bootstrap::Int = 100)

An pairwise independence test based on the [`JointDistanceDistribution`](@ref)
(Amigó & Hirata, 2018)[^Amigo2018].

## Description

In Amigó & Hirata (2018), a one-sided t-test is used to determine whether the
joint distance distribution is positively biased. Due to the number of elements
(`measure.B`) in the distribution being low, we instead use a bootstrapped
estimate of the mean of the joint distance distribution, and compute the p-value
as the proportion of bootstrapped estimates of the mean larger than `measure.μ`.

[^Amigo2018]:
    Amigó, José M., and Yoshito Hirata. "Detecting directional couplings from multivariate
    flows by the joint distance distribution." Chaos: An Interdisciplinary Journal of
    Nonlinear Science 28.7 (2018): 075302.
"""
struct JointDistanceDistributionTest{M <: JointDistanceDistribution, R} <: IndependenceTest
    measure::M
    rng::R
    n_bootstrap::Int

    function JointDistanceDistributionTest(
            measure::M = JointDistanceDistribution();
            rng::R = Random.default_rng(),
            n_bootstrap::Int = 100) where {M, R}
        new{M, R}(measure, rng, n_bootstrap)
    end
end

# Performs a one-sided t-test to see if the joint distance distribution is skewed above
# `measure.μ`, which is the hypothetical mean of the joint distance
# distribution under the null (defaults to `0.0`).

struct JDDTestResult2{V, P}
    Δmeans::V
    pvalue::P
end

function Base.show(io::IO, r::JDDTestResult2)
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
        H₀: "The input variables are independent"
        H₁: "The input variables are dependent"
        ----------------------------------------------------------------------------------
        Bootstrapped Δmeans: $(r.Δmeans)
        p-value:             $(r.pvalue)
          $α005
          $α001
          $α0001\
        """
        )
end

function independence(test::JointDistanceDistributionTest, x, y)
    (; measure, rng, n_bootstrap) = test
    Δjdd = jdd(measure, x, y)
    # B is usually quite low, and data Δjdd not guaranteed to be normal. Therefore,
    # bootstrap the mean and compute a p-value as #(bᵢ) <= μ for i = 1:n, where
    # bᵢ is one of n boostrapped estimates of the mean.
    # Amigó and Hirata use a OneSampleTTest.
    μΔjdd_bootstrapped = bootstrap(mean, x; n = n_bootstrap, rng)
    pval = count(μΔjdd_bootstrapped .<= measure.μ) / length(μΔjdd_bootstrapped)
    return JDDTestResult2(μΔjdd_bootstrapped, pval)
    #return OneSampleTTest(Δjdd, measure.μ)
end
