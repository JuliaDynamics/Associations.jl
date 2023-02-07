using Random: MersenneTwister
rng = MersenneTwister(1234)
x, y = rand(rng, 3000), rand(rng, 3000)
m = JointDistanceDistribution(D = 5, B = 10)
test =  JointDistanceDistributionTest(m)
@test test isa JointDistanceDistributionTest
r = independence(test, x, y)
@test r isa CausalityTools.JDDTestResult

# Don't reject null at significance level (1 - α) when there is no coupling.
α = 0.01
@test pvalue(r) > α 

# Reject null at significance level (1 - α) when there is coupling
z = y .+ x
r = independence(test, z, x)
@test pvalue(r) < α
