using Test
using Associations 
using StableRNGs

rng = StableRNG(123)

function ar3(n::Int, rng = StableRNG(123))
    x = zeros(n)
    y = zeros(n)
    z = zeros(n)
    x[1:3] .= randn(rng, 3)
    y[1:3] .= randn(rng, 3)
    z[1:3] .= randn(rng, 3)

    for i = 4:n
        x[i] = 0.6*sqrt(2)*x[i-1] - 0.71*x[i-2] + randn(rng)
        y[i] = 0.7*x[i-2] + 0.71*y[i-2] + randn(rng)
        z[i] = 0.67*z[i-2] - 0.6*y[i-1] - 0.2*z[i-3] + randn(rng)
    end
    return x, y, z
end

α = 0.05

ηTf = 1
embedding = EmbeddingTE(; dS = 2, dT = 2, dC = 2, ηTf)
measure = TEShannon(; embedding)

# ----------------------------------------------------------------
# Transfer entropy
# ----------------------------------------------------------------
sys = system(Logistic4Chain(; rng))
x, y, z, w = columns(first(trajectory(sys, 50, Ttr = 10000)))
est = Zhu1(TEShannon())
test = LocalPermutationTest(est, nshuffles = 2)
@test independence(test, x, z, y) isa LocalPermutationTestResult

sys = system(Logistic4Chain(; rng))
x, y, z, w = columns(first(trajectory(sys, 50, Ttr = 10000)))
est = EntropyDecomposition(TEShannon(), Kraskov(k = 3))
test = LocalPermutationTest(est, nshuffles = 19)
@test independence(test, x, z, y) isa LocalPermutationTestResult


# For the dedicated estimators, we actually test the outcome on longer timeseries.
# This is because the transfer entropy based local permutation test implemented 
# here doesn't appear in the literature. It is new, so we need to verify that it works.
dedicated_estimators = [Lindner(measure, k=10), Zhu1(measure, k=10)]
@testset "LocalPermutationTest with TEShannon + dedicated TE estimator $(typeof(estimator).name.name)" for estimator in dedicated_estimators
    x, y, z = ar3(200, rng)

    independence_test = LocalPermutationTest(estimator; nshuffles = 19, rng = rng)
    # x and z should be independent given y 
    # (so we shouldn't be able to reject the null, i.e. pvalue >= α)
    @test independence(independence_test, x, z, y).pvalue >= α

    # x and y should NOT be independent given z
    # (so we SHOULD be able to reject the null, i.e. pvalue < α)
    @test independence(independence_test, x, y, z).pvalue < α

    # A test with noise (all variables should be conditionally independent)
    # (so we shouldn't be able to reject the null, i.e. pvalue >= α)
    x, y, z = randn(rng, 100), randn(rng, 100), randn(rng, 100)
    @test independence(independence_test, x, z, y).pvalue >= α
    @test independence(independence_test, x, y, z).pvalue >= α

    # Only conditional analyses are possible with dedicated `TransferEntropyEstimator`s,
    # meaning that we need three inputs. Pairwise analyses don't work.
    @test_throws ArgumentError independence(independence_test, x, y)
end

nondedicated_estimators = [FPVP(), GaussianCMI(), 
    EntropyDecomposition(TEShannon(), Kraskov()), 
    EntropyDecomposition(TEShannon(), PlugIn(Shannon()), CodifyVariables(ValueHistogram(2)))
]
@testset "LocalPermutationTest with TEShannon + non-dedicated estimator $(typeof(estimator).name.name)" for estimator in nondedicated_estimators
    x, y, z = ar3(50, rng)

    independence_test = LocalPermutationTest(estimator; nshuffles = 19, rng = rng)
    @test independence(independence_test, x, z, y) isa LocalPermutationTestResult
end

