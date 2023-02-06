using Random: MersenneTwister
rng = MersenneTwister(1234)
x, y = rand(rng, 1000), rand(rng, 1000)
test =  JointDistanceDistributionTest()
@test test isa JointDistanceDistributionTest
r = independence(test, x, y)
@test r isa CausalityTools.JDDTestResult

# Don't reject null at significance level (1 - α) when there is no coupling.
α = 0.0001
@test pvalue(r) > α 

# Reject null at significance level (1 - α) when there is coupling
x = randn(1000)
y = 0.2randn(1000) .+ x
r = independence(test, x, y)
@test pvalue(r) < 0.0001
