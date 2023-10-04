using Test
using CausalityTools 
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
dedicated_estimators = [Lindner(k=10), Zhu1(k=10)]
@testset "LocalPermutationTest with TEShannon + dedicated TE estimator $estimator" for estimator in dedicated_estimators
    x, y, z = ar3(500, rng)

    independence_test = LocalPermutationTest(measure, estimator; nshuffles = 100, rng = rng)
    # x and z should be independent given y 
    # (so we shouldn't be able to reject the null, i.e. pvalue >= α)
    @test independence(independence_test, x, z, y).pvalue >= α

    # x and y should NOT be independent given z
    # (so we SHOULD be able to reject the null, i.e. pvalue < α)
    @test independence(independence_test, x, y, z).pvalue < α

    # A test with noise (all variables should be conditionally independent)
    # (so we shouldn't be able to reject the null, i.e. pvalue >= α)
    x, y, z = randn(rng, 500), randn(rng, 500), randn(rng, 500)
    @test independence(independence_test, x, z, y).pvalue >= α
    @test independence(independence_test, x, y, z).pvalue >= α

    # Only conditional analyses are possible, meaning that we need three inputs.
    # Pairwise analyses won't work, because only two inputs are given.
    @test_throws ArgumentError independence(independence_test, x, y)
end

nondedicated_estimators = [FPVP(), GaussianMI(), Kraskov(), ValueHistogram(2)]
@testset "LocalPermutationTest with TEShannon + non-dedicated estimator $estimator" for estimator in nondedicated_estimators
    x, y, z = ar3(50, rng)

    independence_test = LocalPermutationTest(measure, estimator; nshuffles = 100, rng = rng)
    @test independence(independence_test, x, z, y) isa LocalPermutationTestResult
end

